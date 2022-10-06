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
  prsRepoMeta <- load_metadata_prsrepo()
  nb <- nrow(prsRepoMeta)

  prsRepoMeta <- prsRepoMeta[!c(prsRepoMeta[["snpRes"]] == snpres),]

  save(prsRepoMeta, file=fs::path(Sys.getenv("PRS_REPO"), "metadata.RDS"))
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






