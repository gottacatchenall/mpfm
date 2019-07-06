load_req_libraries = function(){
    library(ggplot2)
    library(tidyr)
    library(dplyr)
    library(progress)
    library(RColorBrewer)
    library(tools)
    library(subprocess)
    library(future)
    library(gridExtra)
}

# ===================================================
# Merge data frame from each run into single df
# ===================================================
merge_dfs = function(all, this_dir){
    if (!is.null(all)) {
        ret = bind_rows(all, this_dir)
        return(ret)
    }
    return(this_dir)
}

# ===================================================
# Read individual data files
# ===================================================
read_params_file = function(metadata_path, id){
    metadata_table = read.table(metadata_path, sep=",")
    rownames(metadata_table) = metadata_table[,1]
    names = metadata_table[,1]
    vals = metadata_table[,2]

    metadata = data.frame(matrix(ncol=length(names)))
    colnames(metadata) = names
    metadata[1,] = t(vals)
    metadata$id = id
    return(metadata)
}

read_data = function(dir, single=F){
    load_req_libraries()

    # Path names
    indiv_pops_path = "individual_populations.csv"
    pairwise_pops_path = "pairwise_populations.csv"
    pop_dist_path = "individual_populations.csv"
    genome_path = "genome.csv"
    
    params_path = "params.ini"
    diskernel_path = "diskernel.ini"

    # Data frames
    pairwise_pop_data = NULL
    population_data = NULL
    metadata = NULL
    genome_data = NULL
    run_id = 0

    pb = progress_bar$new(total= length(list.files(dir)))

    direcs = NULL
    if (single){
      direcs = c(dir)
      print(direcs)
    }
    else{
      direcs = list.files(dir)
    }
    
    for (direc in direcs){ 
      abs_path = ""
      if (single){
        abs_path = direc
      }
      else{
        abs_path = paste(dir, direc, sep="/")
      }
        pb$tick()

        # Get paths for this dir
        this_dirs_params_path = paste(abs_path, "/", params_path,  sep="")
        this_dirs_individual_pops_path = paste(abs_path, "/", indiv_pops_path, sep="")
        this_dirs_pairwise_pops_path = paste(abs_path, "/", pairwise_pops_path, sep="")
        this_dirs_genome_path = paste(abs_path, "/", genome_path, sep="")
      
        # Read this dir data
        tryCatch({
            print(run_id)
            #this_dirs_pararms = read_metadata_file(this_dirs_metadata_path, run_id)
            this_dirs_indiv_pops = read.csv(this_dirs_individual_pops_path, row.names=NULL, stringsAsFactors=F, na.strings="unknown")
            this_dirs_indiv_pops$id = rep(run_id, nrow(this_dirs_indiv_pops))
            colnames(this_dirs_indiv_pops) = c("generation","pop1","x","y","w_mean","prop_of_k","effective_migration","prop_of_loci_fixed","mean_polymorphism_ct_per_locus","ef0", "blank","sample_size","id")
            this_dirs_indiv_pops[] <- lapply(this_dirs_indiv_pops, function(x) as.numeric(as.character(x)))
            
            this_dirs_pairwise_pop_data = read.csv(this_dirs_pairwise_pops_path, stringsAsFactors=F, na.strings="unknown")
            this_dirs_pairwise_pop_data$id = rep(run_id, nrow(this_dirs_pairwise_pop_data))
            this_dirs_pairwise_pop_data[] <- lapply(this_dirs_pairwise_pop_data, function(x) as.numeric(as.character(x)))
            
            this_dirs_genome_data = read.csv(this_dirs_genome_path, stringsAsFactors=F, na.strings="unknown")
            this_dirs_genome_data$id =  rep(run_id, nrow(this_dirs_genome_data))
            this_dirs_genome_data[] <- lapply(this_dirs_genome_data, function(x) as.numeric(as.character(x)))
            
            this_dirs_metadata = read_params_file(this_dirs_params_path, run_id)
            this_dirs_metadata$id =  rep(run_id, nrow(this_dirs_metadata))
            
            
            
            
            # Merge this dir with the rest
            pairwise_pop_data = merge_dfs(pairwise_pop_data, this_dirs_pairwise_pop_data)
            population_data = merge_dfs(population_data, this_dirs_indiv_pops)
            genome_data = merge_dfs(genome_data, this_dirs_genome_data)
            metadata = merge_dfs(metadata, this_dirs_metadata)
        },
        error = function(e){
            print("errored")
          print(e)
        }
        )
        run_id = run_id + 1
    }

    return(list(pairwise_pop_data, population_data, genome_data, metadata))
}

