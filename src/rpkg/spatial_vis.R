get_population_locations = function(dataframe){
    one_timestep = subset(dataframe, generation==0)
    locations = data.frame(pop=one_timestep$population, x=one_timestep$x, y=one_timestep$y)
    return(locations)
}

get_resistance_values = function(dataframe){
    one_timestep = subset(dataframe, generation==0)
    resistance_values = data.frame(pop1=one_timestep$pop1, pop2=one_timestep$pop2, resistance=one_timestep$resistance_val)
    return(resistance_values)
}

draw_lines_between_points = function(pop_locations, resistance_values, decay_str=1.5, mig_prop=1){
    df = data.frame(matrix(ncol=6))
    
    A = subset(pop_locations, pop == -2)$x
    B = subset(pop_locations, pop == -1)$x
  
    pop_locations = subset(pop_locations, pop >= 0)
    
    x_vals = pop_locations$x
    y_vals = pop_locations$y
    
    n_pops = length(x_vals) - 1
    
    m2 = 1 / (B-A)

    line_ct = 1
    grp_ct = 0 
    for (pt_ct in seq(0, n_pops-1)){
        x1 = subset(pop_locations, pop == pt_ct)$x
        y1 = subset(pop_locations, pop == pt_ct)$y
        
        for (pt_ct2 in seq(pt_ct+1, n_pops-1)){
            if (pt_ct != pt_ct2){
                x2 = subset(pop_locations, pop == pt_ct2)$x
                y2 = subset(pop_locations, pop == pt_ct2)$y
                
                if (length(x) & length(y)){
                    color = 4
                    
                    m1 = (y2 - y1)/(x2-x1)
                    
                    x = (m1 * x1 - m2 * A - y1) / (m1 - m2)
                    y = (x - x1)*m1 + y1
                    
                    # Check conditions on whatever
                    loX = min(x1, x2)
                    hiX = max(x1, x2)
                    loY = min(y1, y2)
                    hiY = max(y1, y2)
                    
                    #print(paste(pt_ct, pt_ct2))
                    #print(paste(loX, hiX, loY, hiY, x,y))
                    
                    resist = subset(resistance_values, pop1 == pt_ct & pop2 == pt_ct2)$resistance
                    #print(paste("p1:", pt_ct, " p2:", pt_ct2))
                    #print(resist)
                    if (length(resist) & resist > 0){
                        #print(paste("does exist", "p1:", pt_ct, " p2: ", pt_ct2))
                        color = 8
                    }
                 
                    dist = sqrt((x2-x1)^2 + (y2-y1)^2)
                    
                    # opac = exp(-1*decay_str*dist)
                    opac = 1.0
                    width = mig_prop
                    
                    df[line_ct,] = c(x1,y1, grp_ct, opac, width, color)
                    line_ct = line_ct + 1
                    df[line_ct,] = c(x2,y2, grp_ct, opac, width, color)
                    line_ct = line_ct + 1
                    grp_ct = grp_ct + 1
                }
            }
        }
    }
    colnames(df) = c("x", "y", "grp", "opac", "width", "type")
    
    df[line_ct,] = c(A, 0, grp_ct, 1.0, 1.5, 2)
    line_ct = line_ct + 1
    df[line_ct,] = c(B, 1, grp_ct, 1.0, 1.5, 2)
    
    plt = ggplot(df, aes(x,y, group=grp)) + geom_line(alpha=df$opac, size=df$width, color= factor(df$type))+ theme(aspect.ratio=1) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line())+ geom_point(size=3.0) + xlim(0,1)
    return(plt)
}

plot_road_fragmentation = function(data){
    pop_data = data[[1]]
    pairwise_data = data[[2]]
    loci_data = data[[3]]
    metadata = data[[4]]
    
    pop_locations =  get_population_locations(pop_data)
    resistance_values = get_resistance_values(pairwise_data)
    
    
    draw_lines_between_points(pop_locations, resistance_values)
    
    # ggplot(df, aes(x,y, group=grp)) + geom_line(alpha=df$opac, size=df$width, color= factor(df$type))+ theme(aspect.ratio=1) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line())+ geom_point(size=1.5) + xlim(0,1)
    #ggplot(pop_locations, aes(x,y group=grp)) + 
}

