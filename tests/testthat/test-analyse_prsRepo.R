create_demo_files <- function(dirpath, string){
  dirfiles <- fs::dir_create(fs::path(dirpath, string))
  files <- fs::path(dirfiles, paste0(string,c(".ma", ".parRes", ".snpRes")))
  fs::file_create(files)
}

create_demo_meta <- function(dir) {
  strings <- c("bip2019","t1d", "scz2022", "an2019", "EA3", "height")
  args <- list(comment = "a gwas",source ="gwascatalogue", source_link="ftp://randomweb.com", population_prevalence = 0.02)

  mdf <- dplyr::bind_rows(
    dplyr::tibble(name = sample(strings, 1),ncase = 100000, ncontrol = 10000, n = 0, other=list(args)),
    dplyr::tibble(name = sample(strings, 1),ncase = 100000, ncontrol = 10000, n = 0, other=list(args)),
    dplyr::tibble(name = sample(strings, 1),ncase = 100000, ncontrol = 10000, n = 0, other=list(args)),
    dplyr::tibble(name = sample(strings, 1),ncase = 100000, ncontrol = 10000, n = 0, other=list(args))
  )
  save(mdf, file = paste0(dir, "/metadata.RDS"))
}





# x <- function(){
#   strings <- c("bip2019","t1d", "scz2022", "an2019", "EA3", "height")
#   # create a temp directory
#   prsrepo <- withr::local_tempdir(tmpdir = "prsRepo")
#   withr::local_envvar("PRS_REPO" = prsrepo)
#   # create lots of fake snpres directories
#   purrr::walk(strings, create_demo_files, dirpath=prsrepo)
#   create_demo_meta(prsrepo)
#   get_pgs_data()
#
#
# }
# dff <- x()
#
#
# dff %>%
#   dplyr::mutate(source = purrr::map(other, "source"),
#                 pop_prev = purrr::map(other, "population_prevalence")
#                 ) %>%
#   dplyr::select(-other) %>%
#   dplyr::mutate(dplyr::across(where(is.list), unlist)) %>%
#   dplyr::mutate(sample_pop = ncase/(ncase+ncontrol)) %>%
#   mutate(neff = prsTools::neff(ncase, ncontrol),)
