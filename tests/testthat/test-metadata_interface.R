

test_that("Errors if directory exists", {
  dir <- withr::local_tempdir(tmpdir = "repo")

  expect_error(setup_repository(dir))
})

test_that("creates a directory", {
  dir <- withr::local_tempdir()
  dir2 <- paste0(dir, "/", "test")
  withr::local_envvar("PRS_REPO"= dir2)
  setup_repository()
  expect_true(fs::dir_exists(dir2))
})

test_that("creates the metadata file", {
  dir <- withr::local_tempdir(pattern = "prsRepo")
  fs::dir_delete(dir)
  withr::local_envvar("PRS_REPO"= dir)
  setup_repository()
  expect_true(fs::file_exists(fs::path(dir, "metadata.RDS")))

})


test_that("remove_all_extensions removes all extensions", {
  fp <- "long/gnarly/file/path/my_gwas.snpres.ma.vcf.gz"

  expect_equal(remove_all_extensions(fp), "my_gwas")
})



