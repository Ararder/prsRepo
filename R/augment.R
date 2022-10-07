ldsc_munge <- function(infile, out){
  glue::glue(
    "$ldsc/munge_sumstats.py ",
    "--sumstats {infile} ",
    "--out {out} ",
    "--merge-alleles $ldsc/w_hm3.snplist",

  )

}
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
subset <- c("bip2021_noUKB", "asd2019")
run_ldsc <- function(subset) {

  tmpdir <- fs::dir_create(Sys.getenv("LD_TEMP"))
  info <- get_snpres_info() %>%
    dplyr::mutate(samp = ncase/(ncase+ncontrol)) %>%
    dplyr::select(snpres, samp, pop = pop_prev)

  inpath <- fs::dir_ls(Sys.getenv("PRS_REPO"), glob ="*ma.gz", recurse = TRUE)
  names <- fs::path_ext_remove(fs::path_ext_remove(fs::path_file(inpath)))

  to_run <-
    dplyr::tibble(snpres=names, inpath) %>%
    dplyr::inner_join(info) %>%
    dplyr::mutate(out = paste0(tmpdir, snpres)) %>%
    dplyr::mutate(run_munge = ldsc_munge(infile = inpath, out = out)) %>%
    dplyr::mutate(run_h2 = ldsc_h2(infile = out,samp = samp, pop = pop))

  # subset if needed.
  if(!missing(subset)){
    to_run <- dplyr::filter(to_run, snpres %in% subset)
  }

  writeLines(c(to_run[["run_munge"]], to_run[["run_h2"]]), "/nfs/home/arvhar/pipelines/gctb/jobs/update_ldsc.sh")
  system("chmod 700 /nfs/home/arvhar/pipelines/gctb/jobs/update_ldsc.sh")
  system("/nfs/home/arvhar/pipelines/gctb/jobs/update_ldsc.sh")
}

update_ldsc <- function() {
  df <- load_metadata_prsrepo()

  # update <- df %>%
  #   dplyr::filter(dplyr::if_all(c(7:15), ~.=="NA"))

  update <- df %>%
    dplyr::filter(dplyr::if_all(c(7:15), is.na))

  run_ldsc(update[["snpRes"]])


}
