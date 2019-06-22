load_network_packages = function(){
    library(ggnet)
    library(igraph)
    library(sna)
    library(network)
}

pairwise_ld_data_to_adj_matrix = function(dataframe, threshold){
    # turn it into an adjacency matrix
    
    matrix_size = length(levels(factor(dataframe$loci1))) + 1
    
    adjacency_matrix = matrix(nrow=matrix_size, ncol=matrix_size)
    nrows = nrow(dataframe)
    
    for (i in seq(1, matrix_size)){
        adjacency_matrix[i,i] = 0
    }
    
    for (i in seq(1, nrows)){
        this_row = dataframe[i,]
        # add 1 because Rs indexing
        x = this_row$loci1 + 1
        y = this_row$loci2 + 1
        if (this_row$ld > threshold){
            adjacency_matrix[x,y] = 1
            adjacency_matrix[y,x] = 1
        }
        else{
            adjacency_matrix[x,y] = 0
            adjacency_matrix[y,x] = 0
        }
    }
    
    return(adjacency_matrix)
}

plot_ld_network_over_time = function(dataframe, ld_threshold, circle=F){
    # subset by gen
    generations = levels(factor(dataframe$generation))
    
    list_of_plots = list()
    
    # need list of all potential interactions 
    
    for (gen in generations){
        this_gen = subset(dataframe, generation == gen)
        this_gen_adj = pairwise_ld_data_to_adj_matrix(this_gen, ld_threshold)
        this_gen_network = network(this_gen_adj)
        if (circle){
            plt = ggnet2(this_gen_network, node.size = 1, mode="circle", color="degree") + labs(title=gen) + theme(legend.position = "none")
        }
        else {
            plt = ggnet2(this_gen_network, node.size = 1, color="degree") + labs(title=gen) + theme(legend.position = "none")
        }
        list_of_plots = c(list_of_plots, list(plt))
    }
    
    #grid.arrange(list(list_of_plots)
    do.call(grid.arrange, c(list_of_plots))
}

plot_ld_network_by_pop = function(dataframe, ld_threshold){
}
