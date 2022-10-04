



# add snpRes
## nametag
### common name
### <pmid>
### ncase
### ncontrol
### n
### date_added





# collect parRes



# run LDSC




# validate phenotypes




# list all available phenotypes



# create a directory at specified location
setup_repository <- function(filepath) {
  if(fs::dir_exists(filepath)) stop("Directory already exists")

  fs::dir_create(filepath, recurse = TRUE)
  fs::file_create(fs::path(filepath, "metadata.tsv"))


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

add_snpres <- function(dirpath, ...) {
  stopifnot(fs::dir_exists(dirpath))
  files <- fs::dir_ls(dirpath)
  required <- c("\\.ma$", "\\.snpRes$", "\\.parRes$")
  for(s in required) stopifnot(any(stringr::str_detect(files, s)))


  if(fs::path_file(dirpath) %in% fs::dir_ls(Sys.getenv("PRS_REPO"))) stop("A folder with that name already exists in PRS_REPO")

  add_metadata_row(dirpath, ...)

  fs::dir_copy(path = dirpath, new_path = Sys.getenv("PRS_REPO"))


}

add_metadata_row <- function(dirpath,...) {
  args(...)
  #dirpath of new folder to add
  name <- fs::path_file(dirpath)

  df <- read_repo_meta()
  new_row <- dplyr::tibble(
    snpRes={{name}},
    common_name = {{ common_name }}, pmid = {{ pmid}}, ncase = {{ ncase }},
    ncontrol = {{ ncontrol }}, n = {{ n}}, date= Sys.Date(), comment = {{ comment }}
    )


  readr::write_tsv(dplyr::add_row(df, new_row), fs::path(Sys.getenv("PRS_REPO"), "metadata.tsv"))


}
read_repo_meta <- function() {
  readr::read_tsv(fs::path(Sys.getenv("PRS_REPO"), "metadata.tsv"))
}



setup_repo_meta <- function() {
  empty <- dplyr::tibble(
    snpRes="",
    common_name = "", pmid ="", ncase = 0,
    ncontrol = 0, n = 0, date= Sys.Date(), comment = "",
    .rows = 0

  )
  readr::write_tsv(empty, fs::path(Sys.getenv("PRS_REPO"), "metadata.tsv"))
}
