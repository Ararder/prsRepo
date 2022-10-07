utils::globalVariables(
  c("V1", "Mean", "n", "pct_nonzero", "ncase", "ncontrol", ".", "snpRes",
    "n_eff", "pop_prev", "h2", "n_snps", "liability_h2_sbayes", "source_link",
    "snpres"))
#' Collects the parRes output from SbayesR
#'
#' @return a tibble
#' @export
#'
#' @examples \dontrun{
#' collect_parres()
#' }
collect_parres <- function() {
  dir <- Sys.getenv("PRS_REPO")
  files <- fs::dir_ls(dir, recurse = TRUE, glob = "*.parRes")


  names <- fs::path_file(fs::path_dir(files))
  suppressWarnings(purrr::map2_df(files, names, get_info))

}


get_info <- function(path, name) {
  parres <- dplyr::tibble(data.table::fread(path))
  h2 <- parres %>%
    dplyr::filter(V1 == "hsq") %>%
    dplyr::pull(Mean)


  n_snps <- parres %>%
    dplyr::filter(stringr::str_detect(V1, "NumSnp")) %>%
    dplyr::summarise(n=sum(Mean)) %>% dplyr::pull(n)
  pct <- parres %>%
    dplyr::filter(stringr::str_detect(V1, "NumSnp")) %>%
    dplyr::mutate(pct =  Mean/sum(Mean)) %>%
    dplyr::filter(V1 != "NumSnp1") %>%
    dplyr::select(bin = V1, n = Mean, pct) %>%
    dplyr::summarise(pct_nonzero = sum(pct)) %>% dplyr::pull(pct_nonzero)

  dplyr::tibble(
    name = {{ name }},
    h2 = h2,
    n_snps = n_snps,
    pct_nonzero = pct
  )
}


#' Calculates liability scale h2 and expected variance explained by PRS
#' Checks all sumstats in PRS_REPO, and extracts key data. WHere possible
#' calculates the expected variance explained in out of sample PRS, and converts
#' the sbayesR heritability estimate to the liability scale.
#' @return a tibble
#' @export
#'
#' @examples \dontrun{
#' get_snpres_info()
#' }
get_snpres_info <- function() {
  mdf <- load_metadata_prsrepo()
  parres <- collect_parres()

  df <- dplyr::mutate(mdf,pop_prev = create_variable(mdf, "population_prevalence")) %>%
    dplyr::mutate(n_eff =  ifelse(n == 0, prsTools::neff(ncase, ncontrol), n)) %>%
    dplyr::inner_join(parres, by = c("snpRes" = "name")) %>%
    dplyr::mutate(sample_prev = ncase/(ncase+ncontrol))

  # catch errors using purrr safely
  safe_lia_h2 <- purrr::safely(prsTools::liability_scale_h2)
  safe_er2 <- purrr::safely(prsTools::e_r2)

  df <- dplyr::mutate(df, liability_h2_sbayes = safe_lia_h2(df[["pop_prev"]], df[["sample_prev"]], df[["h2"]])[["result"]])

  # have to catch cases where function errors (whhich returns a null), and replace with NA
  expected_r2 <- purrr::map2(df[["liability_h2_sbayes"]],df[["n_eff"]], safe_er2, m = 50000) %>%
    purrr::map("result") %>%
    purrr::modify_if(., is.null, ~NA_real_)

  dplyr::mutate(df, expected_r2=unlist(expected_r2)) %>%
    dplyr::select(snpres=snpRes,n_eff, ncase, ncontrol, pop_prev, h2, liability_h2_sbayes, expected_r2,n_snps.x)


}




create_variable <- function(tbl, variable) {
  # variable is a name of list entry. Extract if from the other column
  raw_entries <- purrr::map(tbl[["other"]], {{ variable }})
  # mutate null entries to NA
  null_entries <- purrr::reduce(purrr::map(raw_entries, is.null), c)
  purrr::reduce(purrr::modify_if(raw_entries, null_entries, ~NA_real_), c)
}





#' Extract information about each sumstat used to get a snpRes file
#'
#' @return a tibble
#' @export
#'
#' @examples \dontrun{
#' get_sumstat_info()
#' }
get_sumstat_info <- function() {
  vars <- c("comment","source", "source_link")
  mdf <- load_metadata_prsrepo()
  extracted <- purrr::map(vars, create_variable, tbl = mdf) %>%
    purrr::reduce(dplyr::bind_cols) %>%
    purrr::set_names(vars)
  dplyr::bind_cols(snpres = mdf[["snpRes"]], extracted, date =mdf[["date"]]) %>%
    dplyr::select(snpres, date, source, comment, source_link)


}






#' Extract information from ldsc log file
#' Extracts relevant information from the .log output of LDSC. reads in all *.log
#' files found in dir
#' @param dir filepath for directory
#'
#' @return a tibble
#' @export
#' @examples \dontrun{
#' collect_ldsc("/path/to/ldsc_logs")
#' }
#'
collect_ldsc <- function(dir){


  snpres <- fs::path_ext_remove(fs::path_file(fs::dir_ls(dir, glob = "*.log")))
  safe_munge <- purrr::safely(munge_ldsc_log)

  safe_output <- purrr::map(fs::dir_ls(dir, glob = "*.log"), safe_munge)


  na_tib <- dplyr::tibble(
    lambda_gc = "NA",
    mean_chi2 = "NA",
    liability_h2_ldsc = "NA",
    liability_h2_ldsc_se = "NA",
    n_snps = "NA",
    ratio = "NA",
    ratio_se = "NA",
    intercept = "NA",
    intercept_se = "NA"
  )

  safe_output %>%
    purrr::map("result") %>%
    purrr::modify_if(., is.null, ~na_tib) %>%
    purrr::reduce(dplyr::bind_rows) %>%
    dplyr::bind_cols(snpres=snpres, .)
}










munge_ldsc_log <- function(path){

  lambda_gc <- function(string) {
    lambda_gc <- stringr::str_subset(string, "Lambda GC: ") %>%
      stringr::str_remove("Lambda GC: ")

    mean_chi2 <- stringr::str_subset(string, "Mean Chi") %>%
      stringr::str_remove("Mean Chi^2: ") %>%
      stringr::str_split(" ") %>% .[[1]] %>%
      .[[3]]

    c(lambda_gc, mean_chi2)

  }
  lia <- function(string) {
    res <- stringr::str_subset(string, "Total Liability scale h2: ") %>%
      stringr::str_remove("Total Liability scale h2: ") %>%
      stringr::str_split(" ") %>% .[[1]] %>%
      stringr::str_replace("\\(", "") %>%
      stringr::str_replace("\\)", "")

    c(res[1], res[2])
  }
  intercept <- function(string) {
    res <- stringr::str_subset(string, "Intercept: ") %>%
      stringr::str_remove("Intercept: ") %>%
      stringr::str_split(" ") %>% .[[1]] %>%
      stringr::str_replace("\\(", "") %>%
      stringr::str_replace("\\)", "")

    c(res[1], res[2])

  }
  ratio <- function(string) {


    subset <-  stringr::str_subset(string, "Ratio: ")

    # No numeric ratio given if ratio is > 0.
    if(length(subset) == 0) {
      res <- c()
      subset <- stringr::str_subset(string, "Ratio")
      res[1] <- subset
      res[2] <- subset
    } else {
      res <- subset %>%
        stringr::str_remove("Ratio: ") %>%
        stringr::str_split(" ") %>% .[[1]] %>%
        stringr::str_replace("\\(", "") %>%
        stringr::str_replace("\\)", "")
    }

    c(res[1], res[2])
  }
  n_snps <- function(string) {
    stringr::str_subset(string, "After merging with regression SNP LD") %>%
      stringr::str_remove("After merging with regression SNP LD, ") %>%
      stringr::str_remove(" SNPs remain")

  }

  x <- readLines(path)
  dplyr::tibble(
    lambda_gc = lambda_gc(x)[1],
    mean_chi2 = lambda_gc(x)[2],
    liability_h2_ldsc = lia(x)[1],
    liability_h2_ldsc_se = lia(x)[2],
    n_snps = n_snps(x),
    ratio = ratio(x)[1],
    ratio_se = ratio(x)[2],
    intercept = intercept(x)[1],
    intercept_se = intercept(x)[2]
  )


}



