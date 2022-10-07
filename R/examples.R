# # install_github("ararder/prsRepo")
library(prsRepository)
# # dir <- "/scratch/tmp/arvhar/ldsc_temp"
# #
#
ldsc <- collect_ldsc("/scratch/tmp/arvhar/ld_temp") %>% dplyr::rename(snpRes = snpres)
# get_snpres_info()
# get_sumstat_info()

