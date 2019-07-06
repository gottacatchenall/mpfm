


setwd("~/Projects/mpfmR")
source('./mpfm.R')

# how do we do pure neutral vs not
# N_FITNESS_LOCI_PER_EF = 0
# N_NEUTRAL ADJUST 
param_dict = list("N_INDIVIDUALS"=c(1000, 2000, 5000), "N_FITNESS_LOCI_PER_EF"=c(0, 10, 30))
#run_mpfm(param_dict, mpfm_path="~/mpfmR", num_replicates = 5)
create_run_dirs_and_create_lb_file(fixed_n_loci = 60, mpfm_path="/projects/mica5688/mpfm/bin/mpfm", data_dir_path="/scratch/summit/mica5688")
#load_r_packages()
#data = read_data('./data', single=F)
#pairwise_data = data[[1]]
#pop_data = data[[2]]
#genome = data[[3]]
#metadata = data[[4]]

#ggplot(pairwise_data, aes(generation, mean_f_st , group=interaction(pop1,pop2,sample_size), color=sample_size))  + geom_line(alpha=0.1) + geom_smooth(group=sample_size)  + theme(legend.position = 'none')


#ggplot(pop_data, aes(generation, , group=interaction(pop1, sample_size), color=sample_size)) + geom_point() + geom_line()

#ggplot(genome, aes(locus, mean_global_ld, group=sample_size,color=sample_size )) + geom_line() + geom_point() + facet_wrap(. ~ generation)
