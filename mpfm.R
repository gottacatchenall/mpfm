
source_files = function(){
  MPFM.DIR = "~/Projects/mpfmR"
  setwd(MPFM.DIR)
  #setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  source('./src/rpkg/run_mpfm.R')
  source('./src/rpkg/read_data.R')
  source('./src/rpkg/ld_networks.R')
  source('./src/rpkg/spatial_vis.R')
  load_network_packages()
}


source_files()

## Single Run
data = run_mpfm_single()

## Batch Run
param_dictionary = list(DISPERSAL_DECAY_BURN_IN=c(0.5, 3.0), DISPERSAL_DECAY_FRAGMENTATION=c(7.0, 15.0))
run_batch_mpfm(param_dictionary)

data = read_data('./data')
plot_road_fragmentation(data)

pop_data = data[[1]]
pairwise_data = data[[2]]
loci_data = data[[3]]
metadata = data[[4]]


