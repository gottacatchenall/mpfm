
group_by_params = function(dataframe, metadata, group_by=c("id")){
    num_params_to_group_by = length(group_by)
    param_vals = data.frame(matrix(nrow=1, ncol=num_params_to_group_by))
    colnames(param_vals) = group_by
    
    for (group in group_by){
        levs = levels(factor(metadata[[group]]))
        param_vals[seq(1,length(levs)),group] = levs
    }
    
    param_combos =  expand.grid(param_vals)
    dataframe$grp = rep("", nrow(dataframe))
    for (grp in seq(1, nrow(param_combos))){
        param_row = param_combos[grp,]
        if (!is.null(param_row)){
            this_param_met = metadata
            grp_str = NA
            for (param in group_by){
                this_group_param_val = param_row[[param]]
                if (!is.na(grp_str)){
                    grp_str = paste(grp_str, this_group_param_val, sep=",")
                }
                else{
                    grp_str = this_group_param_val
                }
                this_param_met = subset(this_param_met, this_param_met[[param]] == this_group_param_val)
            }
            
            
            ids = unique(this_param_met$id)
            print(ids)
            for (id_val in ids){
                dataframe[which(dataframe$id == id_val), "grp"] = grp_str
                dataframe[which(dataframe$id == id_val), "grp"] 
            }
        }
    } 
    return(dataframe)
}

plot_by = function(dataframe, x, y){
   ggplot(dataframe, aes(x, y, color=grp, group=interaction(id,grp))) + geom_line() + geom_point() + facet_wrap(. ~ grp)
}

