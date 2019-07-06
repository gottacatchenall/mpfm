

## Go to the right palce 
setwd("~/Projects/mpfmR")
source('./mpfm.R')

## LOCAL
create_run_dirs_and_create_lb_file(param_dict, num_replicates = 5, fixed_n_loci = 60, mpfm_path="./bin/mpfm", data_dir_path="./data", lb_data_dir = "/scratch/summit/mica5688", lb_file_path="lbfile")

# ================================================================
#  stuff to run
#  batch on summit
# ================================================================
param_dict = list("N_INDIVIDUALS"=c(1000, 2000, 5000), "N_FITNESS_LOCI_PER_EF"=c(0, 10, 30), "N_POPULATIONS"=c(10,20,50))
create_run_dirs_and_create_lb_file(param_dict, num_replicates = 50, fixed_n_loci = 60, mpfm_path="/projects/mica5688/mpfm/bin/mpfm", data_dir_path="./data", lb_data_dir = "/scratch/summit/mica5688", lb_file_path="lbfile")




# ================================================================
#  deal w/ data out
# ================================================================
load_r_packages()
data = read_data('~/new_mpfm_summit_data/data', single=F)
pairwise_data = data[[1]]
pop_data = data[[2]]
genome = data[[3]]
metadata = data[[4]]

# crap here we are again doing it in a script 

df = group_by_params(pop_data, metadata, group_by=c("N_INDIVIDUALS", "N_FITNESS_LOCI_PER_EF"))

ggplot(subset(df, sample_size==1.0), aes(generation, mean_polymorphism_ct_per_locus, group=interaction(grp, id, pop1, sample_size), color=factor(pop1))) + geom_point(size=0.1) + geom_line() + facet_wrap(. ~ grp)




#ggplot(pairwise_data, aes(generation, mean_f_st , group=interaction(pop1,pop2,sample_size), color=sample_size))  + geom_line(alpha=0.1) + geom_smooth(group=sample_size)  + theme(legend.position = 'none')


#ggplot(pop_data, aes(generation, , group=interaction(pop1, sample_size), color=sample_size)) + geom_point() + geom_line()

#ggplot(genome, aes(locus, mean_global_ld, group=sample_size,color=sample_size )) + geom_line() + geom_point() + facet_wrap(. ~ generation)
