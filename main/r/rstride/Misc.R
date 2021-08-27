############################################################################# #
#  This file is part of the Stride software. 
#  It is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by 
#  the Free Software Foundation, either version 3 of the License, or any 
#  later version.
#  The software is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  You should have received a copy of the GNU General Public License,
#  along with the software. If not, see <http://www.gnu.org/licenses/>.
#  see http://www.gnu.org/licenses/.
#
#
#  Copyright 2020, Willem L
############################################################################# #
#
# HELP FUNCTIONS FOR rSTRIDE PRE- AND POST-PROCESSING                       
#
############################################################################# #

#.rstride$set_wd()  #DEVELOPMENT: to set the work directory as the latest stride install dir 
#.rstride$load_pd() #DEVELOPMENT: to retrieve the latest project directory (project_dir)

# load required R package
require(doParallel,quietly = TRUE)
require(simid.rtools,quietly = TRUE)

if(!(exists('.rstride'))){
  .rstride <- new.env()
}

# ############################# #
# ## COMMAND LINE MESSAGES   ####
# ############################# #

# terminate rStride
.rstride$cli_abort <- function(message='')
{
  if(nchar(message)>0)
    smd_print(message,WARNING=T)
  
  smd_print('!! TERMINATE rSTRIDE CONTROLLER !!',WARNING=T)
}

# get system memory state
.rstride$print_system_memory_info <- function(par_nodes_info = NA){
  
  if (!is.na("par_nodes_info") && (Sys.getpid() == par_nodes_info$pid_master || 
                                   Sys.getpid() == par_nodes_info$pid_slave1)) {
  
    cmd_macos <- 'top -l 1 -s 0 | grep PhysMem'
    cmd_linux <- 'grep MemFree /proc/meminfo '
    
    if(Sys.info()['sysname'] == 'Darwin'){
      smd_print(system(cmd_macos,intern = T))
    }
    
    if(Sys.info()['sysname'] == 'Linux'){
      # print in kB
      # smd_print(system(cmd_linux,intern = T)) 
      
      # print in GB
      mem_out <- system(cmd_linux,intern = T)
      mem_out <- gsub('MemFree:','',mem_out)
      mem_out <- gsub(' ','',mem_out)
      mem_out <- gsub('kB','',mem_out)
      smd_print('MemFree:',round(as.numeric(mem_out)/1e6),'GB')
    }  
  }
  
}

############################# #
## PROJECT SUMMARY         ####
############################# #

.rstride$load_project_summary <- function(project_dir, stringsAsFactors=FALSE){
  
  # check if project_dir exists
  if(.rstride$dir_not_present(project_dir)){
    stop('PROJECT DIR NOT PRESENT')
  }
  
  # get the filename in the project dir that contains "summary.csv" (full path)
  project_summary_filename <- file.path(project_dir,dir(project_dir,pattern = 'summary.csv'))
  
  # check if the project summary file can be found
  if(length(project_summary_filename) == 0 || !file.exists(project_summary_filename)){
    stop('SUMMARY FILE NOT FOUND')
  }
  
  # read the csv file
  project_summary          <- read.table(project_summary_filename,sep=',',header=T,stringsAsFactors = stringsAsFactors)
  
  # set NA as character
  project_summary[is.na(project_summary)] <- 'NA'
  
  # return the data.frame
  return(project_summary)
}

############################# #
## OPEN PDF STREAM         ####
############################# #

.rstride$create_pdf <- function(project_dir,file_name,width=7,height=7){
  
  # get run_tag
  run_tag           <- basename(project_dir)
  
  # get file name with path
  file_name_path    <- file.path(project_dir,paste0(run_tag,'_',file_name))
  
  # check extension and add if not present
  if(!grepl('.pdf',file_name_path)){
    file_name_path <- paste0(file_name_path,'.pdf')
  }
  
  # open pdf stream
  pdf(file_name_path,width,height)
  
}

############################# #
## OPEN JPG STREAM         ####
############################# #
.rstride$create_jpg <- function(project_dir,file_name,width=7,height=7){
  
  # get run_tag
  run_tag           <- basename(project_dir)
  
  # get file name with path
  file_name_path    <- file.path(project_dir,paste0(run_tag,'_',file_name))
  
  # check extension and add if not present
  if(!grepl('.jpg',file_name_path)){
    file_name_path <- paste0(file_name_path,'.jpg')
  }
  
  # check X11 environment
  if(nchar(Sys.getenv("DISPLAY"))==0){
    options(bitmapType='cairo')
  }
  
  # open pdf stream
  jpeg(file_name_path,width,height,units='in',res=100)
  
  return(file_name_path)
}

############################# #
## CREATE EXPERIMENT TAG   ####
############################# #

# create experiment tag
.rstride$create_exp_tag <- function(i_exp){
  
  # create experiment tag with leading 0's
  exp_tag <- paste0('exp',sprintf("%04s", i_exp))
  
  # solve issue with spaces instead of 0's on linux (leibniz)
  exp_tag <- gsub(' ','0',exp_tag)
  
  # return
  return(exp_tag)
}


############################# #
## XML FUNCTIONS           ####
############################# #

# list_config <- config_disease
# root_name <- 'disease'
# output_prefix <- 'sim_output'
# Save a list in XML format with given root node
.rstride$save_config_xml <- function(list_config,root_name,filename){
  
  # setup XML doc (to add prefix)
  xml_doc = newXMLDoc()
  
  # setup XML root
  root <- newXMLNode(root_name, doc = xml_doc)
  
  # add list info
  smd_listToXML(root, list_config)
  
  # xml prefix
  xml_prefix <- paste0(' This file is part of the Stride software [', format(Sys.time()), ']')
  
  # check filename XML extension, and add if not present
  if(!grepl('\\.xml',filename)){
    paste0(filename,'.xml')
  }
  
  # save as XML,
  # note: if we use an XMLdoc to include prefix, the line break dissapears...
  # fix: http://r.789695.n4.nabble.com/saveXML-prefix-argument-td4678407.html
  cat( saveXML( xml_doc, indent = TRUE, prefix = newXMLCommentNode(xml_prefix)),  file = filename) 
}


# Read XML and reformat numeric values
.rstride$read_config_xml <- function(config_exp_filename,sel_exp=NA){

  # read xml file  
  config_exp <- xmlToList(config_exp_filename)

  # find numeric parameters
  bool_is_numeric <- suppressWarnings(!is.na(as.numeric(config_exp)))
  
  # convert string into numeric
  config_exp[bool_is_numeric] <- as.numeric(config_exp[bool_is_numeric])
  
  if(!is.na(sel_exp)){
    config_exp <- config_exp[sel_exp,]
  }
  
  # return
  return(config_exp)
}


############################# #
## MATRIX OPERATIONS       ####
############################# #
# note: integers are converted to a string with the transpose operation of a mixed matrix...
# so, separate the comparison for numeric and non-numeric types
.rstride$get_equal_rows <- function(f_matrix,f_vector){
  
  # convert NA into character
  f_vector[is.na(f_vector)] <- 'NA'
  f_matrix[is.na(f_matrix)] <- 'NA'
  
  # get numeric columns
  col_numeric      <- unlist(lapply(f_vector,is.numeric))
  
  # compare
  bool_numeric     <- as.logical(colSums(t(f_matrix[,names(f_vector)[col_numeric]]) == c(f_vector[col_numeric])) == sum(col_numeric))
  bool_not_numeric <- as.logical(colSums(t(f_matrix[,names(f_vector)[!col_numeric]]) == c(f_vector[!col_numeric])) == sum(!col_numeric))
  
  # return combined result
  return(bool_numeric & bool_not_numeric)
}

######################################## #
## AGGREGATE EXPERIMENT OUTPUT FILES  ####
######################################## #

.rstride$aggregate_compressed_output <- function(project_dir,get_csv_output){
  
  # load project summary
  project_summary      <- .rstride$load_project_summary(project_dir)
  
  # get output files
  # data_filenames <- unique(dir(file.path(project_summary$output_prefix),pattern='.RData',full.names = T))
  data_filenames <- unique(dir(file.path(project_dir,'exp_all'),pattern='.rds',full.names = T))
  
  # get output types
  data_type_all <- unique(c(names(readRDS(data_filenames[1])),
                            names(readRDS(data_filenames[length(data_filenames)]))))
  
  # loop over the output data types
  data_type <- data_type_all[2]
  for(data_type in data_type_all){

    # # check cluster
    # smd_check_cluster()
     
    # loop over all experiments, rbind
    i_file <- 2
    data_all <- foreach(i_file = 1:length(data_filenames),
                        .combine=.rstride$rbind_fill) %do%
    {
      # get file name
      exp_file_name <- data_filenames[i_file]
      
      # load output data
      data_exp_all    <- readRDS(exp_file_name)
  
      # check if data type present, if not, next experiment
      #if(!data_type %in% names(data_exp_all)){
      if(length(data_exp_all[[data_type]]) == 1 && is.na(data_exp_all[[data_type]])){
        return(NULL)
      }
      
      # select output type
      data_exp        <- data_exp_all[[data_type]]
      
      # for prevalence data, check the number of days
      if(grepl('prevalence',data_type)){

        # create full-size data frame to include the maximum number of days
        data_tmp        <- data.frame(matrix(NA,ncol=max(project_summary$num_days)+2)) # +1 for day 0 and +1 for exp_id
        names(data_tmp) <-  c(paste0('day',0:max(project_summary$num_days)),
                                'exp_id')

        # insert the experiment data
        data_tmp[names(data_exp)] <- data_exp

          # replace the experiment data by the newly constructed data.frame
          data_exp <- data_tmp
        }
        
        # return
        data_exp
      } # end i_file loop
      
    # continue if data_all is not NULL
    if(any(!is.null(data_all))){

      # make data.frame #TODO: contine with data.table 
      data_all <- as.data.frame(data_all)

      # make id's unique => by adding a exp_id tag with leading zero's
      names_id_columns  <- names(data_all)[grepl('id',names(data_all)) & names(data_all) != 'exp_id']
      num_exp_id_digits <- nchar(max(data_all$exp_id))+1
      
      if(length(names_id_columns)>0) {
        for(i_id_column in names_id_columns){
          row_is_id  <- !is.na(data_all[,i_id_column])
          if(i_id_column %in% c('household_id','school_id','college_id','workplace_id','household_cluster_id','collectivity_id')){
            row_is_id <- row_is_id & data_all[,i_id_column] != 0
          } 
          data_all[row_is_id,i_id_column] <- as.numeric(sprintf(paste0('%d%0',num_exp_id_digits,'d'),
                                                                as.numeric(data_all[row_is_id,i_id_column]),
                                                                data_all$exp_id[row_is_id]))
        }
      }
      
      # save
      run_tag <- unique(project_summary$run_tag)
      save(data_all,file=file.path(project_dir,paste0(run_tag,'_',data_type,'.RData')))
      
      if(get_csv_output){
        write.table(data_all,file=file.path(project_dir,paste0(run_tag,'_',data_type,'.csv')),sep=',',row.names=F)
      }
      
      } # end if data_all is not NULL
  } # end data-type loop
}

.rstride$load_aggregated_output <- function(project_dir,file_type,exp_id_opt = NA){
  
  # load project summary
  project_summary <- .rstride$load_project_summary(project_dir)
  
  # get ouput filenames
  dir_files       <- dir(project_dir,full.names = TRUE)
  output_filename <- dir_files[grepl(file_type,dir_files)]
  
  # if the file does not exists, return NA
  if(length(output_filename)==0){
    return(NA)
  }
  
  # load output
  param_name          <- load(output_filename)
  
  # rename parameter
  data_out <- get(param_name)
  
  # selection?
  if(!any(is.na(exp_id_opt))){
    data_out <- data_out[data_out$exp_id %in% exp_id_opt,]
  }
  
  # return
  return(data_out)
  
}


############################### #
## DEFENSIVE PROGRAMMING     ####
############################### #

.rstride$no_return_value <- function(){
  return(invisible())
}

.rstride$dir_not_present <- function(path_dir){
  
  # if directory does not exists, return TRUE (+warning)
  if(!file.exists(paste(path_dir))){
    smd_print('[',paste(match.call(),collapse=' '),'] DIRECTORY NOT PRESENT:',path_dir,WARNING=T)
    return(TRUE)
  } 
  
  # else, return FALSE
  return(FALSE)
}

# check file presence
.rstride$data_files_exist <- function(design_of_experiment = exp_design){
  
  # #TODO: scan automaticaly for .csv or .json or .xml files and check presence
  # c(design_of_experiment)[grepl('\\.csv',c(design_of_experiment))]
  
  # get the unique file names
  file_names <- unique(c(design_of_experiment$age_contact_matrix_file,
                         design_of_experiment$disease_config_file,
                         design_of_experiment$holidays_file,
                         design_of_experiment$population_file))
  
  # add the path to the data folder
  data_dir <- './data'
  file_names <- file.path(data_dir,file_names)
  
  # check the existance of the files
  file_not_exist_bool   <- !file.exists(file_names)
  
  # if any file missing => return FALSE
  if(any(file_not_exist_bool)){
    smd_print('DATA FILE(S) MISSING:', paste(file_names[file_not_exist_bool],collapse = ' '),WARNING=T)
    return(FALSE)
  }  
  
  # else => return TRUE
  return(TRUE)
}

# log level
# check file presence
.rstride$log_levels_exist <- function(design_of_experiment = exp_design){
  
  valid_levels <- design_of_experiment$event_log_level %in% 
    c('None','Incidence','Transmissions','All','ContactTracing')
  
  if(any(!valid_levels)){
    smd_print('INVALID LOG LEVEL(S):', 
              paste(design_of_experiment$event_log_level[!valid_levels],collapse = ' '),
              WARNING=T)
    return(FALSE)
  }  
  
  # else => return TRUE
  return(TRUE)
  
}

# R0
.rstride$valid_r0_values <- function(design_of_experiment = exp_design){
  
  if(any(!is.null(design_of_experiment$r0)))
  {
    
    r0_max <- max(design_of_experiment$r0)
    
    for(disease_config_file in unique(design_of_experiment$disease_config_file)){
      
      # load disease config file
      config_disease    <- xmlToList(file.path('data',disease_config_file))
      
      # get R0 limit    
      fit_r0_limit <- as.numeric(config_disease$label$fit_r0_limit)
      
      # check
      if(r0_max > fit_r0_limit){
        smd_print('INVALID R0 CONFIG VALUE(S):', paste(design_of_experiment$r0,collapse = ' '),paste0('(R0 LIMIT = ',fit_r0_limit,')') ,WARNING=T)
        return(FALSE)
      } 
    } # end for-loop
  }
  
  # else => return TRUE
  return(TRUE)
  
}

# immunity
.rstride$valid_immunity_profiles <- function(design_of_experiment = exp_design){
  
  immunity_profiles <- unique(c(design_of_experiment$immunity_profile,design_of_experiment$vaccine_profile))
  
  # get immunity profile names
  disease_immunity_profiles <- c('None','Random','AgeDependent','Teachers','Cocoon')
  
  # check if given profile names are valid
  if(!all(immunity_profiles %in% disease_immunity_profiles)){
    smd_print('INVALID IMMUNITY PROFILE(S):', paste(immunity_profiles,collapse = ' '),WARNING=T)
    return(FALSE)
  } # end if-compare
  
  # else => return TRUE
  return(TRUE)
  
}

# immunity
.rstride$valid_seed_infected <- function(design_of_experiment = exp_design){
  
  # select unique combinations of population file and seeding rate
  unique_exp_design <- data.frame(population_file= unique(design_of_experiment[,c('population_file')]))
  
  # add the path to the data folder
  data_dir <- './data'
  unique_exp_design$population_file_full <- file.path(data_dir,unique_exp_design$population_file)
  
  # count lines
  unique_exp_design$population_size <- NA
  for(i_file in 1:nrow(unique_exp_design)){
    file_connnection                         <- file(unique_exp_design$population_file_full[i_file]) 
    unique_exp_design$population_size[i_file] <- length(readLines(file_connnection)) - 1              # -1 for header  
    close(file_connnection)
  }
  
  # merge population size with design of experiment parameters
  design_of_experiment <- merge(design_of_experiment,unique_exp_design)
  
  # compare infected seeds with population size... and print warning if needed
  if(any(design_of_experiment$num_infected_seeds > design_of_experiment$population_size)){
    flag_issue <- unique_exp_design$num_infected_seeds > design_of_experiment$population_size
    smd_print('INIALLY INFECTED > POPULATION SIZE:', paste(design_of_experiment[flag_issue,1:2], collapse = ' & initially infected '),WARNING=T)
    return(FALSE)
  }
  
  # else => return TRUE
  return(TRUE)
  
}

.rstride$valid_cnt_param <- function(design_of_experiment = exp_design){
 
  names(exp_design)[grepl('cnt',names(exp_design))]
  
  # contact reduction should be between 0 and 1
  colname_cnt_reduction <- grepl('cnt_reduction',names(exp_design)) & !grepl('cutoff',names(exp_design))
  if(any(exp_design[,colname_cnt_reduction]>1 | exp_design[,colname_cnt_reduction]<0)){
    smd_print('CONTACT REDUCTIONS SHOULD ALL BE BETWEEN [0,1]',WARNING=T)
    return(false)
  }
     
  # compliance delay should be an integer
  colname_compliance <- grepl('compliance_delay',names(exp_design))
  if(any(exp_design[,colname_compliance] != round(exp_design[,colname_compliance]))){
    smd_print('CONTACT COMPLIANCE DELAY IS NOT AN INTEGER',WARNING=T)
    return(FALSE)
  }
  
  if(any(exp_design$cnt_intensity_householdCluster>1 | exp_design$cnt_intensity_householdCluster<0)){
    smd_print('CONTACT INTENSITY IN HOUSEHOLD CLUSTER SHOULD BE BETWEEN [0,1]',WARNING=T)
    return(false)
  }
  
  return(TRUE)
}

.rstride$check_population_contact_combination <- function(exp_design){
  
  bool_population <- any(grepl("collectivity",exp_design$population_file))
  bool_contacts   <- any(grepl("collectivity",exp_design$age_contact_matrix_file))
  
  if((bool_population || bool_contacts) &&
     !(bool_population && bool_contacts)){
   
    smd_print(c("WARNING: no uniformity on collectivities in population and contact input files.\n\t\t\t- ",
                exp_design$population_file,'\n\t\t\t- ',exp_design$age_contact_matrix_file),WARNING = T,FORCED = T)
    
  }
}

############################### #
## DEVELOPMENT FUNCTIONS     ####
############################### #

# load last project_dir
.rstride$load_pd <- function(){
  
  # set most recent build as work directory
  .rstride$set_wd()
  
  # default output dir
  output_dir              <- 'sim_output'
  
  # load directory content (non recursive)
  sim_dirs <- dir(output_dir) 
  
  # exclude .csv files
  sim_dirs <- sim_dirs[!grepl('.csv',sim_dirs)]
  
  # create project_dir (global)
  project_dir <<- file.path(output_dir,sim_dirs[length(sim_dirs)])
  
  # terminal message
  smd_print('SET PROJECT DIR TO ', project_dir)
  
  # load rStride source code
  source('bin/rstride/rStride.R')
  smd_print('rSTRIDE LOADED')
  
}

# set most recent stride install directory as work directory 
#.rstride$set_wd()
.rstride$set_wd <- function(){
  
  stride_dir_tag  <- 'stride-'
  
  # default install directory
  install_dir              <- system('echo $HOME/opt',intern=T)
  
  # if directory does not exists OR contains a 'stride' folder ==>> try VSC SCRATCH directory
  if(!dir.exists(install_dir) && !any(dir(install_dir,pattern = stride_dir_tag))){
    install_dir              <- system('echo $VSC_SCRATCH',intern=T)
  }
  
  # if directory does not exists, abort
  if(!dir.exists(install_dir)){
    
    smd_print('LATEST STRIDE INSTALLATION COULD NOT BE FOUND')
  
  } else {
    # load directory content (non recursive)
    stride_dir_tag  <- 'stride-'
    stride_dirs     <- dir(install_dir,pattern = stride_dir_tag)
    stride_dirs_num <- suppressWarnings(as.numeric(sub(stride_dir_tag,'',stride_dirs)))
    
    # select last directory
    last_stride_dir <- stride_dirs[order(stride_dirs_num,decreasing = T)[1]]
    
    # set work directory
    setwd(file.path(install_dir,last_stride_dir))
    
    # terminal message
    smd_print('NEW WORK DIRECTORY ',file.path(install_dir,last_stride_dir))
  }
  
  
}

.rstride$is_ua_cluster <- function(){
  return(any(grepl('vsc',dir('~',full.names=T))))
}


# cumulative sum, ignoring NA's
cumsum_na <- function(x){
  cumsum(replace_na(x,0))
}

# help function to combine unequal vectors in foreach loop
.rstride$rbind_fill <- function(x,y){
  if(is.null(y)){
    return(x)
  }
  if(length(setdiff(names(x),names(y)))==0){
    return(rbind(x,y))
  } else{
    return(rbind(x,y,fill=TRUE))
  }
}

############################### #
## DESIGN OF EXPERIMENT      ####
############################### #

.rstride$get_full_grid_exp_design <- function(exp_param_list,num_seeds){
  
  # add sequence with all rng seeds
  exp_param_list$rng_seed = seq(num_seeds)
  
  # generate grid with all combinations
  exp_design <- expand.grid(exp_param_list,
                            stringsAsFactors = F)
  
  # add a unique seed for each run
  exp_design$rng_seed <- sample(nrow(exp_design))
  dim(exp_design)
  
  # return
  return(exp_design)
}
  
.rstride$get_lhs_exp_design <- function(exp_param_list,num_experiments,num_rng_seeds=1){
  
  # get number of values per parameter
  num_param_values <- unlist(lapply(exp_param_list,length))
  
  # select the parameters with at least 2 values
  sel_param     <- names(num_param_values[num_param_values>1])
  
  # set rng seed for LHS 
  rng_seed_lhs <- round(mean(exp_param_list$num_infected_seeds))
  set.seed(ifelse(is.na(rng_seed_lhs),1234,rng_seed_lhs))
    
  # setup latin hypercube design (and add parameter names)
  lhs_design <- data.frame(randomLHS(num_experiments,length(sel_param)))
  names(lhs_design) <- sel_param
  
  # rescale LHS to given parameter range
  for(i_param in sel_param){
    param_range <- range(exp_param_list[i_param])
    lhs_design[,i_param] <- lhs_design[,i_param] * diff(param_range) + param_range[1]
  }
  
  # fix for parameters that are a (non-decimal) number
  sel_param_num <- names(lhs_design)[grepl('num',names(lhs_design)) | grepl('compliance_delay',names(lhs_design))]
  lhs_design[,sel_param_num] <- round(lhs_design[,sel_param_num])
  
  # copy lhs design into 'exp_design' and add other parameters
  exp_design <- lhs_design
  exp_param_names <- names(exp_param_list)
  exp_param_names <- exp_param_names[!exp_param_names %in% sel_param]
  i_param = exp_param_names[1]
  for(i_param in exp_param_names){
    exp_design[,i_param] <- exp_param_list[i_param]
  }
  
  # copy lhs design if num_rng_seeds > 1
  if(num_rng_seeds>1){
    row_ind <- rep(1:nrow(exp_design),each=num_rng_seeds)
    exp_design <- exp_design[row_ind,]
    dim(exp_design)
  }
  
  # add a unique seed for each run
  exp_design$rng_seed <- sample(nrow(exp_design))
  dim(exp_design)
  
  # return
  return(exp_design)
}



