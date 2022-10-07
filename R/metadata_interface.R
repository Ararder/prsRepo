utils::globalVariables("prsRepoMeta")


# create a directory at specified location
setup_repository <- function() {
  if(Sys.getenv("PRS_REPO") == "") stop("No environment variable for PRS_REPO")
  filepath <- Sys.getenv("PRS_REPO")
  if(fs::dir_exists(filepath)) stop("Directory already exists")

  fs::dir_create(filepath, recurse = TRUE)

  prsRepoMeta <-
    dplyr::tibble(
      snpRes="ts",ncase = 0, ncontrol = 0, n = 0, date= Sys.Date(), other = list(list(t="x", y="t"))
    )

  save(prsRepoMeta, file=fs::path(Sys.getenv("PRS_REPO"), "metadata.RDS"))


}

#' Loads the info on all snpRes files in PrsRep
#'
#' @return a tibble
#' @export
#'
#' @examples \dontrun{
#' load_metadata_prsrepo()
#' }
load_metadata_prsrepo <- function(){
  load(paste0(Sys.getenv("PRS_REPO"), "/metadata.RDS"))
  return(dplyr::tibble(prsRepoMeta))
}

# remove all extensions

remove_all_extensions <- function(path) {
  name <- fs::path_file(path)
  while(stringr::str_detect(name, "\\.")){
    name <- fs::path_ext_remove(name)
  }
  name

}

# add a sumstat folder to repository

#' Add a snpRes folder to the prsRepo file
#'
#' @param dirpath filepath for the folder to copy over
#' @param ncase n cases
#' @param ncontrol n controls
#' @param n total n (for continous traits)
#' @param ... Other arguments
#'
#' @return nothing
#' @export
#'
#' @examples \dontrun{
#' add_snpres("~/snpres/gwas", common_name = "mdd gwas", pmid=1231313,
#' ncase=14000, control = 14000, comment = "first mdd GWAS by PGC")}
add_snpres <- function(dirpath, ncase=0, ncontrol=0, n=0, ...) {
  stopifnot(fs::dir_exists(dirpath))
  files <- fs::dir_ls(dirpath)
  required <- c("\\.ma$", "\\.snpRes$", "\\.parRes$")
  for(s in required) stopifnot(any(stringr::str_detect(files, s)))


  if(fs::path_file(dirpath) %in% fs::dir_ls(Sys.getenv("PRS_REPO"))) stop("A folder with that name already exists in PRS_REPO")

  add_metadata_row(dirpath, ncase = {{ ncase }}, ncontrol = {{ ncontrol }}, n = {{ n}}, ...)

  fs::dir_copy(path = dirpath, new_path = Sys.getenv("PRS_REPO"))

  # gzip ma and snpres file inside copied folder
  new_path <- fs::path(Sys.getenv("PRS_REPO"), fs::path_file(dirpath))
  snpres <- fs::dir_ls(new_path, glob ="*.snpRes")
  ma     <- fs::dir_ls(new_path, glob ="*.ma")
  system(paste("gzip", snpres))
  system(paste("gzip", ma))



}

add_metadata_row <- function(dirpath, ncase=0, ncontrol=0, n=0, ...) {
  #dirpath of new folder to add
  name <- fs::path_file(dirpath)
  args <- list(...)

  load(paste0(Sys.getenv("PRS_REPO"), "/metadata.RDS"))
  prsRepoMeta <- dplyr::add_row(prsRepoMeta, snpRes={{name}}, ncase = {{ ncase }}, ncontrol = {{ ncontrol }}, n = {{ n}}, date= Sys.Date(), other = list(args))
  save(prsRepoMeta, file=fs::path(Sys.getenv("PRS_REPO"), "metadata.RDS"))

}

#' Remove a row from the metadata file
#'
#' @param snpres name of the snpres file that you want to delete (folder name in
#' prsRepo)
#'
#' @return nothing
#' @export
#'
#' @examples \dontrun{
#' delete_metadata_row("scz2022")
#' }
#'
delete_metadata_row <- function(snpres) {
  if(is.integer(snpres)) {
    prsRepoMeta <- load_metadata_prsrepo()
    prsRepoMeta <- prsRepoMeta[-snpres,]
    save(prsRepoMeta, file=fs::path(Sys.getenv("PRS_REPO"), "metadata.RDS"))
  } else {
    stopifnot(is.character(snpres))
    prsRepoMeta <- load_metadata_prsrepo()
    prsRepoMeta <- prsRepoMeta[!c(prsRepoMeta[["snpRes"]] == snpres),]
    save(prsRepoMeta, file=fs::path(Sys.getenv("PRS_REPO"), "metadata.RDS"))
  }

}


#' Update an existing row in the metadata
#'
#' @param snpres name of snpres file
#' @param new_row the new updated row
#'
#' @return nothing
#' @export
#'
#' @examples \dontrun{
#'
#' update_metadata_row("scz2022", new_row)
#' }
update_metadata_row <- function(snpres, new_row) {
  prsRepoMeta <- load_metadata_prsrepo()

  prsRepoMeta[c(prsRepoMeta[["snpRes"]] == snpres),] <- new_row

  save(prsRepoMeta, file=fs::path(Sys.getenv("PRS_REPO"), "metadata.RDS"))
}

utils::globalVariables(
  c("V1", "Mean", "n", "pct_nonzero", "ncase", "ncontrol", ".", "snpRes",
    "n_eff", "pop_prev", "h2", "n_snps", "liability_h2_sbayes", "source_link",
    "snpres"))



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





