load_requirements = function(){
    require(plotly)
}

mk_mpfm_instance_directory = function(data_dir_path, instance_dir_name){
    # Creates data_dir_path if it doesn't exist
    dir.create(data_dir_path, showWarnings = FALSE)
    
    if (is.null(instance_dir_name)){
        now = Sys.time()
        now_str = format(now,'%B_%d_%Y__%H:%M:%S')
        instance_dir_name =  now_str
    }
    path = paste(data_dir_path, instance_dir_name, sep="/")
    dir.create(path)
    return(path)
}

create_instance_param_file = function(params, instance_dir){
    param_file_path = paste(instance_dir, "params.ini", sep="/")
    params_to_write = subset(params, value != "")
    write.table(params_to_write, file= param_file_path, sep=",", row.names = F, col.names= F, quote=F)
}


run_mpfm_instance = function(params, data_dir_path="./data", delete_instance_dir_when_done=T){
    instance_dir_name = params$DATA_DIRECTORY
    instance_dir = mk_mpfm_instance_directory(data_dir_path, instance_dir_name)
    
    mpfm_exe_path = "./bin/mpfm"
    
    create_instance_param_file(params, instance_dir)
  
    ## Execute mpfm 
    cmd = paste(mpfm_exe_path, instance_dir)
    print(cmd)
    system(cmd)
    
    # Get data back in correct form
    
    # Aggregate all runs if more than one 
    
    
    # Clean uo run directory 
    if (delete_instance_dir_when_done){
        unlink(instance_dir, recursive=T) 
    }
}


get_parameters = function(){
    param_table = read.csv("./src/param_table.csv")
    param_df = data.frame(matrix(ncol=2, nrow=nrow(param_table)))
    colnames(param_df) = c("parameter", "value")
    param_df$parameter = param_table$Parameter
    param_df$value = param_table$Default
    return(param_df)
}

params_test = get_parameters()

run_mpfm_instance(params_test, delete_instance_dir_when_done = F)

