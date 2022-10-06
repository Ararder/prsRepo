sbayes_call <-  prsTools::construct_template(
  program = "/nfs/home/arvhar/pipelines/gctb/gctb",
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

#' Title
#'
#' @param sumstat filepath for .ma file to be rescaled
#' @param name give the output a name
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
