
source_files = function(){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  source('./src/rpkg/run_mpfm.R')
  source('./src/rpkg/read_data.R')
  source('./src/rpkg/ld_networks.R')
  source('./src/rpkg/spatial_vis.R')
  load_network_packages()
}

run_mpfm_single = function(){
  load_requirements()
  params = get_parameters()
  run_mpfm_instance(params, delete_instance_dir_when_done = F)
  
  this_dir = dirname(rstudioapi::getSourceEditorContext()$path)
  data_directory = paste(this_dir, "data", sep="/")
  data = read_data(data_directory)
  return(data)
}

run_batch_mpfm = function(){
  load_requirements()
}

data = run_mpfm_single()


source_files()
data = read_data('./data')

plot_road_fragmentation(data)

pop_data = data[[1]]
pairwise_data = data[[2]]
loci_data = data[[3]]
metadata = data[[4]]


ggplot(subset(loci_data, pop1 == 0 & pop2 == 3), aes(generation, clustering, group=threshold, color=as.factor(threshold))) + geom_point() + geom_line()

ggplot(subset(loci_data, pop1 == 0 & pop2 == 1), aes(generation, modularity, group=threshold, color=as.factor(threshold))) + geom_point() + geom_line()

# look at dist and ibr vs delta mod and delta clus during fragmentation 

# mod 
this_tmp_df  = loci_data %>% right_join(pairwise_data, by=c("pop1"="pop1","pop2"="pop2", "generation"="generation", "id"="id"))
plt1 = ggplot(subset(this_tmp_df, threshold==0.5 & resistance_val == 0), aes(generation, modularity, color=dist, group=interaction(threshold,pop1, pop2))) + geom_point() + geom_line()
plt2 = ggplot(subset(this_tmp_df, threshold==0.5 & resistance_val == 1), aes(generation, modularity, color=dist, group=interaction(threshold,pop1, pop2))) + geom_point() + geom_line()
grid.arrange(plt1, plt2)

# clus
plt1 = ggplot(subset(this_tmp_df, threshold==0.5 & resistance_val == 0), aes(generation, clustering, color=dist, group=interaction(threshold,pop1, pop2))) + geom_point() + geom_line()
plt2 = ggplot(subset(this_tmp_df, threshold==0.5 & resistance_val == 1), aes(generation, clustering, color=dist, group=interaction(threshold,pop1, pop2))) + geom_point() + geom_line()
grid.arrange(plt1, plt2)

## LD Networks
plot_ld_network_over_time(subset(loci_data, pop1 == 0 & pop2 == 3), 0.1)




### Delta LD corrected for distance
deltas = data.frame(matrix(nrow = 1, ncol=5))
colnames(deltas) = c("resist", "dist", "delta_cl", "delta_mod", "id")
i = 1
for (idv in levels(as.factor(pairwise_data$id))){
  post_frag = subset(this_tmp_df, id == idv & generation == 2280 )
  pre_frag = subset(this_tmp_df, id == idv & generation == 1900 )
  
  delta_cl = post_frag$modularity - pre_frag$modularity
  delta_mod = post_frag$clustering - pre_frag$clustering
  nrows = length(delta_cl)
  if (nrow(post_frag) > 0 && nrow(pre_frag) > 0 && nrows > 0){
    deltas = rbind(deltas, data.frame(resist=post_frag$resistance_val,dist=post_frag$dist, delta_cl=delta_cl, delta_mod=delta_mod, id=rep(idv, nrows)))
    i = i + nrows
  }
}

mod = ggplot(deltas, aes(x=resist, y=delta_mod, group=resist, color=dist)) + geom_violin() + geom_point() 
cl = ggplot(deltas, aes(x=resist, y=(delta_cl), group=resist, color=dist)) + geom_violin(alpha=0.5)+ geom_point()

grid.arrange(mod, cl)


## LD Clustering coeff over time 

# for each pair of pops ij
cutoffs = c(0.01, 0.05, 0.1, 0.25, 0.4)


get_clustering_coeff_over_time(loci_data, 3, 16, cutoffs)
get_clustering_coeff_over_time(loci_data, 8, 16, cutoffs)


get_clustering_coeff_over_time = function(loci_data, p1, p2, cutoffs){
  tmp_data = subset(loci_data, pop1 == p1 & pop2 == p2)
  # make a df to hold data 
  this_pair_of_pops_clustering = data.frame(matrix(ncol=4, nrow=1))
  colnames(this_pair_of_pops_clustering) = c("clus", "mod", "gen", "cutoff") 
  generations = levels(factor(tmp_data$generation))
  
  line = 1
  
  for (gen in generations){
    this_data = subset(tmp_data, generation == gen)
    for (cutoff in cutoffs){
      adj_matrix = pairwise_ld_data_to_adj_matrix(this_data, cutoff)
      g = graph_from_adjacency_matrix(adj_matrix, mode="undirected")
      clustering_coefficient = transitivity(g, isolates = "zero")
      
      mod = modularity(cluster_fast_greedy(g))
      
      this_pair_of_pops_clustering[line,] = c(clustering_coefficient, mod, gen, cutoff)
      
      #else{ 
       # this_pair_of_pops_clustering[line,] = c(0, gen, cutoff)
      #}
      line = line + 1
    }
    # plot time vs. clus for multiple values of cutoff on same plot 
  }
  
  this_pair_of_pops_clustering$clus = as.numeric(this_pair_of_pops_clustering$clus)
  this_pair_of_pops_clustering$cutoff = as.numeric(this_pair_of_pops_clustering$cutoff)
  this_pair_of_pops_clustering$gen = as.numeric(this_pair_of_pops_clustering$gen)
  this_pair_of_pops_clustering$mod = as.numeric(this_pair_of_pops_clustering$mod)
  
  ggplot(this_pair_of_pops_clustering, aes(gen, clus, group=cutoff, color=factor(cutoff))) + geom_line() + geom_point() + labs(title = paste("IBR ; IBD"))

  
  ggplot(this_pair_of_pops_clustering, aes(gen, mod, group=cutoff, color=factor(cutoff))) + geom_line() + geom_point() + labs(title = paste("IBR ; IBD"))
  
}


## LD along genome 
whoops = aggregate(ld ~ generation*loci1, FUN=mean, data=loci_data)
ggplot(whoops, aes(loci1, ld)) + geom_line() + facet_wrap(. ~ generation) 

ggplot(loci_data, aes(generation, ld, color=interaction(loci1, loci2), group=interaction(loci1, loci2))) + geom_point() + geom_line() + theme(legend.position="none") 

### Linkge and FST comparison, change during frag

ggplot(subset(pairwise_data), aes(generation, f_st, color=dist, group=interaction(pop1, pop2, id))) + geom_point(alpha=0.3, size=0.2) + geom_line(alpha=0.5, size=1.9) + scale_color_distiller(palette="YlOrRd") + theme_bw() + scale_x_continuous(breaks = seq(0, 2400, by = 600), limits = c(0,2400), expand = c(0, 0)) + scale_y_continuous(limits = c(0,0.5), expand = c(0, 0))+ theme(text = element_text(size=16, family="Menlo"), aspect.ratio = 1, axis.ticks.length=unit(.45, "cm"))  + xlab("Generation") + ylab("F_ST") + labs(title="Pairwise F_ST") 

ggplot(subset(pairwise_data), aes(generation, ld, color=dist, group=interaction(pop1, pop2, id))) + geom_point(alpha=0.7) + geom_line(alpha=0.7, size=1.9)  + scale_color_distiller(palette="YlOrRd") + theme_bw() + scale_x_continuous(limits = c(0,2300), expand = c(0, 0)) + scale_y_continuous(limits = c(0,0.3), expand = c(0, 0))+ theme(text = element_text(size=20, family="Menlo"), aspect.ratio = 1)+ xlab("Generation") + ylab("LD") + labs(title="Pairwise Linkage Disequilibrium")



### Delta LD corrected for distance
deltas = data.frame(matrix(nrow = 1, ncol=5))
colnames(deltas) = c("resist", "dist", "delta_fst", "deltald", "id")
i = 1
for (idv in levels(as.factor(pairwise_data$id))){
  post_frag = subset(pairwise_data, id == idv & generation == 2300 )
  pre_frag = subset(pairwise_data, id == idv & generation == 1900 )
  
  delta_ld = post_frag$ld - pre_frag$ld 
  delta_fst = post_frag$f_st - pre_frag$f_st
  nrows = length(delta_fst)
  if (nrow(post_frag) > 0 && nrow(pre_frag) > 0 && nrows > 0){
    deltas = rbind(deltas, data.frame(resist=post_frag$resistance_val,dist=post_frag$dist, delta_fst=delta_fst, deltald=delta_ld, id=rep(idv, nrows)))
    i = i + nrows
  }
}


fst_model = lm(delta_fst ~ resist * dist , data = deltas)
ld_model = lm(deltald ~ resist * dist , data = deltas)

ggplot(deltas, aes(resist, dist, z=deltald)) + geom_point(aes(color=deltald))


f_st = ggplot(deltas, aes(x=resist, y=delta_fst, group=resist, color=dist)) + geom_violin() + geom_point() + ylim(-0.03, 0.1)
ld = ggplot(deltas, aes(x=resist, y=(deltald), group=resist, color=dist)) + geom_violin(alpha=0.5)+ geom_point()+ ylim(-0.03, 0.1)

grid.arrange(f_st, ld)
