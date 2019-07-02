

setwd("~/Projects/mpfmR")
source('./mpfm.R')


param_dict = list(DISPERSAL_DECAY_BURN_IN=c(3.0,0.5), DISPERSAL_DECAY_FRAGMENTATION=c(5.0, 10.0))
run_mpfm(param_dict, num_replicates = 5)
