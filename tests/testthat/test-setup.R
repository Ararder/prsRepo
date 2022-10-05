
setup_repo <- function()
  withr::local_envvar()
# if filepath already exists, setup_repository does not work.
test_that("Errors if directory exists", {
  dir <- withr::local_tempdir(tmpdir = "repo")

  expect_error(setup_repository(dir))
})

test_that("creates a directory", {
  dir <- withr::local_tempdir()
  dir2 <- paste0(dir, "/", "test")
  setup_repository(dir2)
  expect_true(fs::dir_exists(dir2))
})

test_that("creates the metadata file", {
  dir <- withr::local_tempdir()
  setup_repository(fs::path(dir, "new"))
  expect_true(fs::file_exists(fs::path(dir, "new/metadata.tsv")))

})


test_that("remove_all_extensions removes all extensions", {
  fp <- "long/gnarly/file/path/my_gwas.snpres.ma.vcf.gz"

  expect_equal(remove_all_extensions(fp), "my_gwas")
})

test_that("add_snpres reads in the metadata file", {

  dir <- withr::local_tempdir()
  file <- fs::path(dir, "metadata.tsv")
  fs::file_create(file)
  withr::local_envvar("PRS_REPO" = dir)
  fs::file_exists(file)

  expect_equal(read_repo_meta(), readr::read_tsv(file))

})

