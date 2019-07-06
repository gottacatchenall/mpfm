
param_dict = list(SAMPLE_SIZE=c(-1, 0.1, 0.2, 0.5), N_INDIVIDUALS=c(1000, 2000, 5000), RANDOM_SEED=c(5))
#run_mpfm(param_dict, num_replicates = 5)

run_mpfm = function(param_dict, num_replicates = 1, populations=NULL, random_pops_each_run=F, mpfm_path = "~/Projects/mpfm",IBR=F){
    if (!housekeeping_checks(mpfm_path=mpfm_path)){
        return;
    }

    treatments = get_treatment_parameter_dfs(param_dict)
    # How to tell if using IBR or IBD?
    if (IBR){
      dis_kernel = setup_IBR()
    }

    # Otherwise, load custom IBR distances

    run_dirs = create_run_directories_for_treatments(treatments, num_replicates, populations=populations)

    proc_handles = init_processes(run_dirs, parallel=F)

}

init_processes = function(run_dirs, parallel=F){
    mpfm_exe = paste(normalizePath(dirname(".")), "bin/mpfm", sep="/")

    if (parallel){
        # do the multithreading
    }

    ## One at a time
    else {
        for (this_dir in run_dirs){
            system(paste(mpfm_exe, this_dir))
            #handle = spawn_process(mpfm_exe, this_dir)
            #process_wait(handle)
        }
    }
}

log_status = function(msg, stdout=T){
  now = Sys.time()
  now_str = format(now,'%Y-%m-%d--%H:%M:%S:')
  if (stdout){
    print(paste(now_str, "  ",msg))
  }
}

create_random_pops = function(){

}

create_run_directories_for_treatments = function(treatments, num_replicates, populations=NULL, same_pops_for_each_treatment=F, data_dir_path = 'data/'){
    paths = c()

    # check the length of treatments to see if its a single run or not
    treat_ct = 0
    for (treatment in treatments){
      if (same_pops_for_each_treatment){
        n_pops = get_parameter_value(treatment, "N_POPULATIONS")
        n_indivs = get_parameter_value(treatment, "N_INDIVIDUALS")
        populations = create_random_populations(n_pops, n_indivs)
      }
        for (rep_ct in seq(1, num_replicates)){
            run_dir_path = paste("treatment", treat_ct, "_rep", rep_ct, sep="")
            full_dir_path = create_run_directory(treatment, data_dir_path, populations=populations, instance_dir_name = run_dir_path)
            paths = c(paths, full_dir_path)
        }
        treat_ct = treat_ct + 1
    }
    return(paths)
}

create_run_directory = function(treatment, data_dir_path, populations=NULL, instance_dir_name = NULL){
    # make sure ./data exists and create it if not?
    #dir.create(data_dir_path, warning=F)


    # this will never be true atm
    if (is.null(instance_dir_name)){
        now = Sys.time()
        now_str = format(now,'%B_%d_%Y__%H:%M:%S')
        instance_dir_name =  now_str
    }

    this_instance_path = paste(data_dir_path, instance_dir_name, sep="")

    this_dir = normalizePath(dirname("."))
    full_path = paste(this_dir, this_instance_path, sep="/")
    dir.create(full_path)

    # if populations is not null, make sure that the df is correct format and write pops.ini. otherwise generate random pops and write pops ini. check that n_pops is the same as pops df
    if (!is.null(populations)){
      create_pops_ini_file(populations, full_path)
      nrow(populations)
      # make sure n_pops in treatment df is same as nrow of pop df
    }
    else{
      n_pops = get_parameter_value(treatment, "N_POPULATIONS")
      n_indivs = get_parameter_value(treatment, "N_INDIVIDUALS")
      pops = create_random_populations(n_pops, n_indivs)
      create_pops_ini_file(pops, full_path)
    }

    bi_dist_decay_str = get_parameter_value(treatment, "DISPERSAL_DECAY_BURN_IN")
    frag_dist_decay_str = get_parameter_value(treatment, "DISPERSAL_DECAY_FRAGMENTATION")



    bi_diskern = get_ibd_diskern(pops, bi_dist_decay_str)
    frag_diskern = get_ibd_diskern(pops, frag_dist_decay_str)

    genome = get_default_genome(treatment)

    create_genome_ini_file(genome, full_path)
    create_genome_ini_file(genome, full_path)
    create_diskern_ini_file(bi_diskern, full_path,  file_name="diskern_pre.ini")
    create_diskern_ini_file(frag_diskern, full_path,  file_name="diskern_post.ini")
    create_params_ini_file(treatment, full_path)
    return(full_path)
}

get_default_genome = function(treatment){
  n_neutral = get_parameter_value(treatment, "N_NEUTRAL_LOCI")
  n_ef = get_parameter_value(treatment, "EF_NUMBER")
  n_loci_per_ef = get_parameter_value(treatment, "N_FITNESS_LOCI_PER_EF")
  n_loci = n_neutral + (n_ef * n_loci_per_ef)
  n_chromo = get_parameter_value(treatment, "N_CHROMOSOMES")

  genome_length = get_parameter_value(treatment, "GENOME_LENGTH")
  init_poly_mean = get_parameter_value(treatment, "INIT_NUM_ALLELES_MEAN")

  map_dist_step = genome_length / n_loci

  genome = data.frame(matrix(ncol=6,nrow=1))
  #genome = na.omit(genome)
  # map distance is distance from 0 on that chromosome
  colnames(genome) = c("locus","selection_str","ef","chromosome","map_distance","init_polymorphism")

  loci_per_ch = ceiling(n_loci / n_chromo)

  chromo_num = 0
  map_ct = 0.0
  for (l in seq(0,n_loci-1)){
    if (l %% loci_per_ch == 0){
      map_ct = 0.0
      chromo_num = chromo_num + 1
    }
    map_ct = map_ct + map_dist_step
    poly_ct = rpois(1, init_poly_mean)

    genome[l,] = c(l,0,-1,chromo_num,map_ct,poly_ct)
  }


  for (ef in seq(0, n_ef-1)){
    this_ef_loci = sample(seq(1,n_loci), n_loci_per_ef)
    print(this_ef_loci)
    for (l in this_ef_loci){
      sel = runif(1)
      genome[l,2] = sel
      genome[l,3] = ef
    }
  }
  return(na.omit(genome))
}

get_ibd_diskern = function(pops, str){
  n_pops = nrow(pops)

  kern = matrix(nrow=n_pops, ncol=n_pops)

  for (i in seq(1, n_pops)){
    i_x = pops[i,1]
    i_y = pops[i,2]

    row_sum = 0

    for (j in seq(1, n_pops)){
      if (i != j){
        j_x = pops[j,1]
        j_y = pops[j,2]
        dist = sqrt((j_y - i_y)^2 + (j_x - i_x)^2)
        kern_val = exp(-1*str*dist)
        kern[i,j] = kern_val
        row_sum = row_sum + kern_val
      }
      if (i == j){
        kern[i,j] = 0
      }
    }

    for (j in seq(1, n_pops)){
      #print('not normalizing dis kern')
       kern[i,j] = kern[i,j] / row_sum
    }
  }

  return(kern)
}

create_random_populations = function(n_pops, n_indivs){
  k = n_indivs / n_pops

  populations = data.frame(matrix(ncol=4))
  colnames(populations) = c("x", "y", "k", "ef0")
  for (i in seq(1, n_pops)){
    x = runif(1)
    y = runif(1)
    populations[i,] = c(x,y,k,1)
  }
  return(populations)
}

create_genome_ini_file = function(genome, full_path){
  path = paste(full_path, "genome.ini", sep="/")
  write.table(genome, file=path, sep=",",row.names=F, col.names=F, quote=F)
}

create_diskern_ini_file = function(diskern, full_path, file_name="diskern.ini"){
  n_pops = nrow(diskern)
  line_ct = 1
  diskern_list = data.frame(matrix(ncol=3))
  colnames(diskern_list) = c("i","j","val")
  for (i in seq(1, n_pops)){
    for (j in seq(1, n_pops)){
      line_val = diskern[i,j]
      diskern_list[line_ct,] = c(i,j,line_val)
      line_ct = line_ct + 1
    }
  }
  path = paste(full_path, file_name, sep="/")
  write.table(diskern_list, file=path, sep=",",row.names=F, col.names=F, quote=F)
}

create_pops_ini_file = function(pops, full_path){
  path = paste(full_path, "pops.ini", sep="/")
  write.table(pops, file=path, sep=",", col.names= F, quote=F)
}

create_params_ini_file = function(treatment, full_path){
  params_to_write = subset(treatment, value != "")

  path = paste(full_path, "params.ini", sep="/")
  write.table(params_to_write, file=path, sep=",", row.names = F, col.names= F, quote=F)
}

create_genomes_ini_file = function(genome, full_path){
  path = paste(full_path, "genomes.ini", sep="/")
  write.table(genome, file=path, sep=",", row.names = F, col.names= F, quote=F)
}

get_treatment_parameter_dfs = function(param_dict, fixed_n_loci = NULL){
    cartesian_product_of_params = expand.grid(param_dict)
    base_df = read_default_params()

    treatment_dfs = list()
    for (i in seq(1,nrow(cartesian_product_of_params))){
      this_set = cartesian_product_of_params[i,]
      this_treatment_param_df = base_df
      for (param in names(this_set)){
        val = this_set[[param]]
        this_treatment_param_df = set_parameter_value(this_treatment_param_df, param, val)
      }
      
      if (!is.null(fixed_n_loci)){
          N_EF = get_parameter_value(this_treatment_param_df, "EF_NUMBER")
          N_LOCI_PER_EF = get_parameter_value(this_treatment_param_df, "N_FITNESS_LOCI_PER_EF")
          n_fitness = N_EF * N_LOCI_PER_EF
          n_neutral = fixed_n_loci - n_fitness
          this_treatment_param_df = set_parameter_value(this_treatment_param_df, "N_NEUTRAL_LOCI", n_neutral)
      }
      
      
      treatment_dfs[[i]] = this_treatment_param_df
    }
    return(treatment_dfs)
}

set_parameter_value = function(param_df, param_name, value){
  index = which(param_df$parameter == param_name)
  if (index){
    param_df[index,2] = value
  }
  return(param_df)
}

get_parameter_value = function(param_df, param_name){
  index = which(param_df$parameter == param_name)
  if (index){
    return(param_df[index,2])
  }
}

read_default_params = function(){
    param_table = read.csv("./src/param_table.csv")
    param_df = data.frame(matrix(ncol=2, nrow=nrow(param_table)))
    colnames(param_df) = c("parameter", "value")
    param_df$parameter = param_table$Parameter
    param_df$value = (param_table$Default)
    return(param_df)
}

create_run_dirs_and_create_lb_file = function(num_replicates = 5, fixed_n_loci=60, mpfm_path="./bin/mpfm", lb_file_path="./lb_file.txt"){
  if (!housekeeping_checks()){
    return;
  }
  treatments = get_treatment_parameter_dfs(param_dict, fixed_n_loci = fixed_n_loci)
  run_dirs = create_run_directories_for_treatments(treatments, num_replicates, populations=NULL, same_pops_for_each_treatment = T)
  lb_file = create_lb_file(run_dirs, mpfm_path, lb_file_path)
}

create_lb_file = function(run_dirs, mpfm_path, lb_file_path){
  for (dir in run_dirs){
    exe_string = paste(mpfm_path, " ", dir, ";", sep="" )
    write(exe_string,file=lb_file_path,append=TRUE)
  }
}



housekeeping_checks = function(mpfm_path="~/Projects/mpfmR"){
    setwd(mpfm_path)
    package_st = load_r_packages()
    compilation_st = compile_mpfmcore()

    if (!package_st || !compilation_st){
        return(0)
    }
    return(1)
}

load_r_packages = function(){
    library(ggplot2)
    library(tidyr)
    library(dplyr)
    library(progress)
    library(RColorBrewer)
    library(tools)
    library(subprocess)
    library(future)
    library(gridExtra)
  
    source('./src/rpkg/read_data.R')
    source('./src/rpkg/spatial_vis.R')
    source('./src/rpkg/read_data.R')

    return(1)
}

compile_mpfmcore = function(){
    # rpkg "serve" has a way to run make files, also remotely???
    # unclear what the best to do here is
    #system('make')
    #a = system('pwd')
    return(1)
}
