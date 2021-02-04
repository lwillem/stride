#!/usr/bin/env Rscript
############################################################################ #
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
#  Copyright 2019, Willem L, Kuylen E & Broeckhove J
############################################################################ #
# 
# Helpfunction(s) to parse the log file(s)
#
############################################################################ #

################################# #
## PARSE LOGFILE               ####
################################# #

"DEVELOPMENT CODE"
if(0==1){
  
  #f_exp_dir <- file.path(output_dir,output_exp_dirs[i_exp])
  f_exp_dir <- file.path(project_dir,'exp0002')
  event_logfile <- file.path(f_exp_dir,'event_log.txt')
  exp_id <- 2;bool_parse_tracing=TRUE
  xx <- parse_event_logfile(event_logfile,2)
  event_logfile <- event_log_filename
  exp_id <- i_exp;bool_parse_tracing=TRUE
}
parse_event_logfile <- function(event_logfile,exp_id,
                                bool_parse_tracing=TRUE)  # reducted transmission output
{

  # terminal message
  cat("PARSING LOGFILE:",event_logfile,fill=TRUE)

  # get LOG categories, by reading file line by line and select 1st column
  data_log_cat <- fread(event_logfile, sep=' ', fill=TRUE,
                           select = 1)
  data_log_cat <- unlist(unique(data_log_cat))

  # initialise output variables
  rstride_out <- list()
  
  # Parse log file using the following tags tags: 
  # - PART    participant info
  # - PRIM    seed infection
  # - TRAN    transmission event
  # - TRAN_M  transmission event (minimal info)
  # - CONT    contact event
  # - VACC    additional immunization
  # - TRACE   contact tracing
  # - UNITEST] universal test log
  # - UNITEST-ISOLATE universal test strategy isolations
  # - NCOM    non-compliance to social distancing measures

  ###################### #
  ## PARTICIPANT DATA ####
  ###################### #
  header_part         <- c('local_id', 'part_age', 'household_id', 'school_id', 'college_id','workplace_id', 
                           'household_cluster_id','collectivity_id',
                           'is_susceptible','is_infected','is_infectious','is_recovered','is_immune',
                           'start_infectiousness','start_symptomatic','end_infectiousness','end_symptomatic',
                           'household_size','school_size','college_size','workplace_size','primarycommunity_size','secundarycommunity_size','is_teleworking')

  rstride_out$data_participants  <- reformat_log_data(event_logfile = event_logfile,
                                                      data_log_cat  = data_log_cat,
                                                      log_cat       = "PART",
                                                      colnames_all  = header_part,
                                                      exp_id        = exp_id)

  
  ####################### #
  ## TRANSMISSION DATA ####
  ####################### #
  if("[TRAN_M]" %in% data_log_cat){
    header_transm       <- c('part_age',
                             'sim_day',
                             'start_infectiousness',
                             'start_symptoms',
                             'end_symptoms')
  } else {
    header_transm       <- c('local_id', 'infector_id','part_age',
                             'infector_age','pool_type','sim_day','id_index_case',
                             'start_infectiousness','end_infectiousness','start_symptoms','end_symptoms',
                             'infector_is_symptomatic','part_rel_infectiousness','part_rel_susceptibility')
  }
  
  rstride_out$data_transmission  <- reformat_log_data(event_logfile = event_logfile,
                                                      data_log_cat  = data_log_cat,
                                                      log_cat       = c("PRIM","TRAN","TRAN_M"),
                                                      colnames_all  = header_transm,
                                                      exp_id        = exp_id)
  
  ###################### #
  ## CONTACT DATA     ####
  ###################### # 
  header_cnt          <- c('local_id', 'part_age', 'cnt_age', 'cnt_home', 'cnt_school', 
                           'cnt_college','cnt_work', 'cnt_prim_comm', 'cnt_sec_comm', 'cnt_hh_cluster', 'cnt_collectivity',
                           'sim_day', 'cnt_prob', 'trm_prob','part_sympt','cnt_sympt')

  rstride_out$data_contacts <- reformat_log_data(event_logfile = event_logfile,
                                                 data_log_cat  = data_log_cat,
                                                 log_cat       = "CONT",
                                                 colnames_all  = header_cnt,
                                                 exp_id        = exp_id)

  
  ###################### #
  ## VACCINATION DATA ####
  ###################### # 
  header_vac          <- c('local_id', 'part_age', 'pool_type', 'pool_id', 'pool_has_infant', 'sim_day')

  rstride_out$data_vaccination <- reformat_log_data(event_logfile  = event_logfile,
                                                    data_log_cat  = data_log_cat,
                                                    log_cat        = "VACC",
                                                    colnames_all   = header_vac,
                                                    exp_id         = exp_id)

  ########################################## #
  ## CONTACT TRACING DATA ####
  ########################################## # 
  header_trace           <- c('local_id', 'part_age', 'is_infected', 'is_symptomatic','pool_type', 
                              'case_id','case_age','sim_day','num_unique_contacts','num_contacts_tested')

  if(bool_parse_tracing){
    rstride_out$data_tracing <- reformat_log_data(event_logfile = event_logfile,
                                                  data_log_cat  = data_log_cat,
                                                  log_cat       = "TRACE",
                                                  colnames_all  = header_trace,
                                                  exp_id        = exp_id)
  }

  
  ########################################## #
  ## UNIVERSAL TESTING: LOG  ####
  ########################################## # 
  header_testing          <- c('sim_day', 'unitest_day_in_sweep')
  
  rstride_out$data_unitest <- reformat_log_data(event_logfile = event_logfile,
                                                data_log_cat  = data_log_cat,
                                                log_cat       = "UNITEST",
                                                colnames_all  = header_testing,
                                                exp_id        = exp_id)
  ########################################## #
  ## UNIVERSAL TESTING: ISOLATION  ####
  ########################################## # 
  header_testing_iso      <- c('pcr_pool_id', 'household_id', 'local_id', 'is_infected', 'isolation_delay', 'sim_day')
  
  rstride_out$data_unitest_iso <- reformat_log_data(event_logfile = event_logfile,
                                                    data_log_cat  = data_log_cat,
                                                    log_cat       = "UNITEST-ISOLATE",
                                                    colnames_all  = header_testing_iso,
                                                    exp_id        = exp_id)
  

  ########################################## #
  ## NON-COMPLIANCE                       ####
  ########################################## # 
  header_ncompliance   <- c('local_id', 'part_age', 
                           'household_id', 'noncomplier_household',
                           'school_id', 'noncomplier_school',
                           'college_id', 'noncomplier_college',
                           'workplace_id',  'noncomplier_workplace',
                           #'hh_cluster_id', 'noncomplier_hh_cluster',
                           'prim_comm_id','noncomplier_prim_comm',
                           'sec_comm_id', 'noncomplier_sec_comm')
  

  rstride_out$ncompliance <- reformat_log_data(event_logfile = event_logfile,
                                               data_log_cat  = data_log_cat,
                                               log_cat       = "NCOM",
                                               colnames_all  = header_ncompliance,
                                               exp_id        = exp_id)
  
  # print CLI message and return
  cat("LOG PARSING COMPLETE",fill=TRUE)
  return(rstride_out)
}

################################# #
## REFORMAT LOG DATA           ####
################################# #
# function to select and reformat log output into a data.table with numeric, booleans and text data
 log_cat      = c("PRIM","TRAN")
# log_cat      = c("'\\[PRIM]'") 
# log_cat      = "PRIM" 
# 
# # log_cat       = "CONT"
# # colnames_all <- header_cnt
log_cat       = "UNITEST-ISOLATE"
# colnames_all  = header_testing_iso
#  colnames_all <- header_transm
# exp_id <- 2
reformat_log_data <- function(event_logfile,data_log_cat,log_cat,colnames_all,exp_id) {

  # adapt log category to ...]
  log_cat <- paste0(log_cat,']')
  
  # check if the given log category is present
  #FIX for unitest: gsub
  if(!any(sapply(log_cat,grepl,data_log_cat))){
    return(NA)
  }
  
  # create grep command based on the given log category/categories
  cmd_grep <- paste(c("grep",log_cat),collapse=' -e ')
  
  # set columns to drop and select final colnames
  colnames_doc    <- c(NA,colnames_all)
  columns_select  <- which(!is.na(colnames_doc))
  colnames_select <- colnames_doc[!is.na(colnames_doc)]
  
  # read file line by line and select the requested lines
  data_log_subset <- fread(cmd=paste(cmd_grep,event_logfile), sep=' ',
                           select = columns_select,
                           col.names = colnames_select)

  # check
  dim(data_log_subset)
  object.size(data_log_subset) / 1e6
  
  # get columns with numeric and boolean values
  colnames_char           <- colnames_select[grepl('type',colnames_select)]
  colnames_boolean        <- colnames_select[grepl('is_',colnames_select)]
  colnames_numeric        <- colnames_select[!colnames_select %in% c(colnames_char,colnames_boolean)]
  
  # specify help functions for 'lapply'
  set_NAs <- function(x){x[x==-1] <- NA; x}
  is_char_true <- function(x){x=='true'}
  remove_tag <- function(x){gsub('.*=','',x)}
 
  # make sure the columns do not contain a parse-tag 
  #TODO: fix in C++ for UNITTEST ?
  if(any(grepl('=',data_log_subset[1, ])))
    data_log_subset[, c(colnames_numeric) := lapply(.SD,remove_tag), .SDcols = colnames_numeric]
  
  # set -1 to NA
  data_log_subset[, c(colnames_select) := lapply(.SD, set_NAs), .SDcols = colnames_select]
  
  # make sure that numeric values are stored as integers
  data_log_subset[, c(colnames_numeric) := lapply(.SD, as.numeric), .SDcols = colnames_numeric]
  
  # convert character booleans 'true' and 'false' into R booleans
  if(length(colnames_boolean)>0)
  data_log_subset[, c(colnames_boolean) := lapply(.SD, is_char_true), .SDcols = colnames_boolean]
  
  # add exp_id
  data_log_subset[,exp_id := exp_id]
  
  # return
  return(data_log_subset)
}


