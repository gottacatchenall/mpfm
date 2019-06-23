
run_mpfm_single = function(){
  load_requirements()
  params = get_parameters()
  dir_path = initialize_mpfm_instance(params)
  
  log_status(paste("Executing single instance of MPFM in", dir_path))
  execute_mpfm_instance(dir_path)
  log_status(paste("MPFM instance in", dir_path, "completed."))
  
  
  this_dir = dirname(rstudioapi::getSourceEditorContext()$path)
  data_directory = paste(this_dir, "data", sep="/")
  data = read_data(data_directory)
  return(data)
}

log_status = function(msg, stdout=T){
  now = Sys.time()
  now_str = format(now,'%B_%d_%Y__%H:%M:%S')
  
  if (stdout){
    print(paste(now_str, "---", msg))
  }
}


# Take in a dictionary of parameters, do each pairwise combo n_rep time
# Param dictionary should be of the form list(PARAM_NAME_1=c(val1, val2,...), PARAM_NAME_2=c(val1, val2))... etc

param_dictionary = list(DISPERSAL_DECAY_BURN_IN=c(0.5, 3.0), DISPERSAL_DECAY_FRAGMENTATION=c(7.0, 15.0))

run_batch_mpfm = function(param_dictionary, n_reps=1){
  load_requirements()
  param_df = get_parameters()
  params_to_set = names(param_dictionary)

  cartesian_product_of_params = expand.grid(param_dictionary)
  treat_ct = 0
  for (i in seq(1,nrow(cartesian_product_of_params))){
    this_set = cartesian_product_of_params[i,]
    this_treatment_param_df = param_df
    for (param in names(this_set)){
      val = this_set[[param]]
      print(val)
      this_treatment_param_df = set_parameter_value(this_treatment_param_df, param, val)
    }
    
    print(this_treatment_param_df)
    for (rep_ct in seq(1,n_reps)){
      this_path = paste("treatment", treat_ct, "_rep", rep_ct, sep="")
      instance_dir_path = initialize_mpfm_instance(this_treatment_param_df, this_path)
    }
    treat_ct = treat_ct + 1
  }
}

# =================================
# Helper Functions
# =================================
set_parameter_value = function(param_df, param_name, value){
  print(param_name)
  print(value)
  index = which(param_df$parameter == param_name)
  if (index){
    param_df[index,2] = value
  }
  print(param_df)
  return(param_df)
}


load_requirements = function(){
    require(plotly)
}

mk_mpfm_instance_directory = function(data_dir_path, instance_dir_name=NULL){
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


initialize_mpfm_instance = function(params, instance_dir_name=NULL, data_dir_path="./data"){
    desc = "Metapopulation Fragmentation Model (mpfm)
    v1.0

    A individual-based model of metapopopulation dynamics in spatiotemporally stochastic environments.

    University of Colorado at Boulder
    Dept. Ecology and Evolutionary Biology | Flaxman Lab
    Michael Catchen | michael.catchen@colorado.edu"


    instance_dir = mk_mpfm_instance_directory(data_dir_path, instance_dir_name)
    
    print(instance_dir)
    
    create_instance_param_file(params, instance_dir)
    return(instance_dir)
}

execute_mpfm_instance = function(instance_dirs, delete_instance_dir_when_done=F){
  for (instance_dir in instance_dirs){
    mpfm_exe_path = "./bin/mpfm"
    cmd = paste(mpfm_exe_path, instance_dir)
    print(cmd)
    system(cmd)
    # Clean uo run directory
    if (delete_instance_dir_when_done){
      unlink(instance_dir, recursive=T)
    }
  }
}


get_parameters = function(){
    param_table = read.csv("./src/param_table.csv")
    param_df = data.frame(matrix(ncol=2, nrow=nrow(param_table)))
    colnames(param_df) = c("parameter", "value")
    param_df$parameter = param_table$Parameter
    param_df$value = (param_table$Default)
    return(param_df)
}

