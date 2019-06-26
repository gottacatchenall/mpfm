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
read_indiv_pop_data = function(indiv_pop_path, demo_path, id){
    demography = read.csv(demo_path)
    indiv_pop = read.csv(indiv_pop_path)

    pop_data_this_run = indiv_pop %>% full_join(demography, by=c("generation", "population"))
    pop_data_this_run[is.na(pop_data_this_run)] = 0
    pop_data_this_run$id = id
    return(pop_data_this_run)
}

read_pairwise_loci_data = function(pairwise_loci_path, id){
    pairwise_by_loci = read.csv(pairwise_loci_path)
    pairwise_by_loci$id = id
    return(pairwise_by_loci)
}

read_pairwise_pop_data = function(pairwise_dist_path, pairwise_ld_path, dispersal_path, id){
    pop_dist = read.csv(pairwise_dist_path)
    pairwise_ld = read.csv(pairwise_ld_path)
    dispersal = read.csv(dispersal_path)

    pairwise_data_this_run = pop_dist %>% right_join(pairwise_ld, by=c("pop1"="population1","pop2"="population2")) %>% left_join(dispersal, by=c("generation", "pop1"="populationA", "pop2"="populationB"))
    pairwise_data_this_run[is.na(pairwise_data_this_run)] = 0
    pairwise_data_this_run$id = id
    return(pairwise_data_this_run)
}


read_metadata_file = function(metadata_path, id){
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

read_data = function(dir){
    load_req_libraries()

    # Path names
    pops_path = "populations.csv"
    pop_dist_path = "pop_distance.csv"
    pairwise_pop_ld_path = "pairwise_ld.csv"
    pairwise_loci_ld_path = "ld_networks.csv"
    dispersal_path = "dispersal.csv"
    demography_path = "demography.csv"
    metadata_path = "params.ini"

    # Data frames
    pairwise_loci_data = NULL
    pairwise_pop_data = NULL
    population_data = NULL
    metadata = NULL
    run_id = 0

    pb = progress_bar$new(total= length(list.files(dir)))

    for (direc in list.files(dir)){
        abs_path = paste(dir, direc, sep="/")
        pb$tick()

        # Get paths for this dir
        this_dirs_metadata_path = paste(abs_path, "/", metadata_path, sep="")
        this_dirs_dispersal_path = paste(abs_path, "/", dispersal_path, sep="")
        this_dirs_pairwise_loci_path = paste(abs_path, "/", pairwise_loci_ld_path, sep="")
        this_dirs_demography_path = paste(abs_path, "/", demography_path, sep="")
        this_dirs_pop_data_path = paste(abs_path, "/", pops_path, sep="")
        this_dirs_pairwise_ld_path = paste(abs_path, "/", pairwise_pop_ld_path, sep="")
        this_dirs_pairwise_dist_path = paste(abs_path, "/", pop_dist_path, sep="")


        # Read this dir data
        tryCatch({
            print(run_id)
            this_dirs_metadata = read_metadata_file(this_dirs_metadata_path, run_id)
            this_dirs_pop_data = read_indiv_pop_data(this_dirs_pop_data_path, this_dirs_demography_path, run_id)
            this_dirs_pairwise_loci_data = read_pairwise_loci_data(this_dirs_pairwise_loci_path, run_id)
            this_dirs_pairwise_pop_data = read_pairwise_pop_data(this_dirs_pairwise_dist_path, this_dirs_pairwise_ld_path, this_dirs_dispersal_path, run_id)

            # Merge this dir with the rest
            pairwise_loci_data = merge_dfs(pairwise_loci_data, this_dirs_pairwise_loci_data)
            pairwise_pop_data = merge_dfs(pairwise_pop_data, this_dirs_pairwise_pop_data)
            population_data = merge_dfs(population_data, this_dirs_pop_data)
            metadata = merge_dfs(metadata, this_dirs_metadata)
        },
        error = function(e){
            print("errored")
        }
        )
        run_id = run_id + 1
    }

    return(list(population_data, pairwise_pop_data, pairwise_loci_data, metadata))
}

get_mig_rate = function(data){
    return(levels(as.factor(subset(pairwise_data, id==id_val)$BASE_MIGRATION_RATE)))
}

get_disp_kern = function(data){
    return(levels(as.factor(subset(pairwise_data, id==id_val)$DISPERSAL_DECAY)))
}
