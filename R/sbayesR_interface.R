utils::globalVariables(
  c("V1", "Mean", "n", "pct_nonzero", "ncase", "ncontrol", ".", "snpRes",
    "n_eff", "pop_prev", "h2", "n_snps", "liability_h2_sbayes", "source_link",
    "snpres"))
# Construct a function that generates a template with the default
# sbayesr arguments
sbayes_call <-  prsTools::construct_template(
  program = Sys.getenv("GCTB"),
  default_args = list(
    sbayes = "R",
    mldm = Sys.getenv("GCTB_LDMATRIX"),
    robust = "",
    pi = "0.95,0.02,0.02,0.01",
    gamma = "0.0,0.01,0.1,1",
    "ambiguous-snp" = "",
    "impute-n" = "",
    "chain-length" = 6000,
    "burn-in" = 1200,
    thin = 10,
    "out-freq" = 10,
    "no-mcmc-bin"= ""
  )
)
sbayes_call_no_impute <-  prsTools::construct_template(
  program = "/nfs/home/arvhar/pipelines/gctb/gctb",
  default_args = list(
    sbayes = "R",
    mldm = Sys.getenv("GCTB_LDMATRIX"),
    robust = "",
    pi = "0.95,0.02,0.02,0.01",
    gamma = "0.0,0.01,0.1,1",
    "ambiguous-snp" = "",
    "chain-length" = 6000,
    "burn-in" = 1200,
    thin = 10,
    "out-freq" = 10,
    "no-mcmc-bin"= ""
  )
)

#' Rescale summary statistics with SbayesR
#' A wrapper and handler for rescaling effect sizes using sbayesR in GCTB.
#' Generate the bash code to run sbayesR, by getting filepaths for environmental
#' variables (GCTB, GCTB_LDMATRIX)
#' @param sumstat filepath for .ma file to be rescaled
#' @param name give the output a name, if not provided, will use filename of
#' filepath, with extension removed.
#' @param impute_n Should sbayesR impute sample size?
#' @param archive send output to Sys.getenv("PRS_ARCHIVE")
#' @param sbatch should the function sbatch the commandline code?
#'
#' @return nothing
#' @export
#'
#' @examples \dontrun{
#' run_sbayes("here/is_my/gwas_sumstat/sumstat.ma", name = "anorexia_2019")
#' }
run_sbayes <- function(sumstat, name, impute_n=TRUE, archive=TRUE, sbatch=FALSE) {
  # if no name is passed, remove extension and use basename for provided filepath
  if(missing(name)) name <- fs::path_ext_remove(fs::path_file(sumstat))
  archive <- Sys.getenv("PRS_ARCHIVE")


  outdir     <- paste0(archive, "/", name)
  snpRes_out <- paste0(outdir,   "/", name)
  slurmfile  <- paste0(outdir,  "/", name,"-%j.out")
  bash_out   <- paste0(outdir,  "/", name, "job.sh")

  fs::dir_create(outdir)
  if(!fs::file_exists(paste0(outdir, "/", fs::path_file(sumstat)))) fs::file_copy(sumstat, outdir)



  header <- c(
    "#!/bin/bash",
    "#SBATCH --time 480",
    "#SBATCH --ntasks 4",
    "#SBATCH --mem 40000",
    paste0("#SBATCH -o ", slurmfile)
  )
  if(impute_n){
    body <- sbayes_call("gwas-summary" = paste0(outdir, "/", fs::path_file(sumstat)), out = snpRes_out)
  } else {
    body <- sbayes_call_no_impute("gwas-summary" = paste0(outdir, "/", fs::path_file(sumstat)), out = snpRes_out)
  }


  writeLines(c(header,body), bash_out)

  if(sbatch) system(paste0("sbatch", " ", bash_out))



}


#' Collects the parRes output from SbayesR
#' Checks all folders recursively starting at the PRS_REPO env variable.
#' Reads in all files matching "*.parRes" --> logfile from sbayesR
#' Extracts relevant data.
#' @return a tibble
#' @export
#'
#' @examples \dontrun{
#' collect_parres()
#' }
collect_parres <- function() {
  files <- fs::dir_ls(Sys.getenv("PRS_REPO"), recurse = TRUE, glob = "*.parRes")
  names <- fs::path_file(fs::path_dir(files))
  suppressWarnings(purrr::map2_df(files, names, get_info))

}


get_info <- function(path, name) {
  # work horse for parsing the parRes file.
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
