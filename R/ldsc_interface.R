utils::globalVariables(c(
  "samp","out","pop"))

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

  # find all the outputs from ldsc
  snpres <- fs::path_ext_remove(fs::path_file(fs::dir_ls(dir, glob = "*.log")))

  # create a safe version for potential error handling of edge cases
  safe_munge <- purrr::safely(munge_ldsc_log)

  # extract all data from ldsc output
  safe_output <- purrr::map(fs::dir_ls(dir, glob = "*.log"), safe_munge)

  # if one output contains null elements, it means something resutled in an
  # error. Add this tibble in that case.
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

  # Extract results, if error, input NA tibble
  safe_output %>%
    purrr::map("result") %>%
    purrr::modify_if(., is.null, ~na_tib) %>%
    purrr::reduce(dplyr::bind_rows) %>%
    dplyr::bind_cols(snpres=snpres, .)
}




munge_ldsc_log <- function(path){
  # A set of helper functions for extracting relevant metrics from LDSC output

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

# helper function for generating a call to ldsc:s munge_sumstats.py
ldsc_munge <- function(infile, out){
  glue::glue(
    "$ldsc/munge_sumstats.py ",
    "--sumstats {infile} ",
    "--out {out} ",
    "--merge-alleles $ldsc/w_hm3.snplist",

  )

}
# helper function for generating a call to ldsc.py --h2
ldsc_h2 <- function(infile, samp, pop) {
  glue::glue(
    "$ldsc/ldsc.py ",
    "--h2 {infile}.sumstats.gz ",
    "--ref-ld-chr $ldsc/eur_w_ld_chr/ ",
    "--w-ld-chr $ldsc/eur_w_ld_chr/ ",
    "--out {infile} ",
    "--samp-prev {samp} ",
    "--pop-prev {pop}"
  )
}

run_ldsc <- function(subset) {
  if(!missing(subset)){
    stopifnot(is.character(subset))
  }
  # output ldsc munged files and output in a temp directory.
  tmpdir <- fs::dir_create(Sys.getenv("LD_TEMP"))

  # get sample prevalence and population prevalence
  info <- get_snpres_info() %>%
    dplyr::mutate(samp = ncase/(ncase+ncontrol)) %>%
    dplyr::select(snpres, samp, pop = pop_prev)

  # list all files ending in *.ma.gz (sumstats in ma format) i prsRepo
  inpath <- fs::dir_ls(Sys.getenv("PRS_REPO"), glob ="*ma.gz", recurse = TRUE)
  # create names from these files
  names <- fs::path_ext_remove(fs::path_ext_remove(fs::path_file(inpath)))

  # bind all these data together, and create calls for munge sumstats and ldsc
  to_run <-
    dplyr::tibble(snpres=names, inpath) %>%
    dplyr::inner_join(info) %>%
    dplyr::mutate(out = paste0(tmpdir, snpres)) %>%
    dplyr::mutate(run_munge = ldsc_munge(infile = inpath, out = out)) %>%
    dplyr::mutate(run_h2 = ldsc_h2(infile = out,samp = samp, pop = pop))

  # if subset variable is passed, only retain those rows matching ids in subset
  if(!missing(subset)){
    to_run <- dplyr::filter(to_run, snpres %in% subset)
  }

  # write a bash script, give permission to run and run it.
  run_script <- paste0(Sys.getenv("TEMP_LD"), "update_ldsc.sh")
  writeLines(c(to_run[["run_munge"]], to_run[["run_h2"]]), run_script)
  system(paste("chmod 700", run_script))
  system(run_script)
}

update_ldsc <- function() {
  # load metadata
  df <- load_metadata_prsrepo()

  # check if any sumstats are missing ldsc rows
  update <- dplyr::filter(df,dplyr::if_all(c(7:15), is.na))

  # run ldsc for missing sumstats
  run_ldsc(update[["snpRes"]])


}



