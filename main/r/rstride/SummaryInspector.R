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
#
#  Copyright 2020, Willem L, Kuylen E & Broeckhove J
############################################################################ #
#
# MODEL SUMMARY EXPLORATION
# - input-output behavior
# - transmission events and context
#
############################################################################ #

inspect_summary <- function(project_dir)
{
  
  # load project output summary
  project_summary    <- .rstride$load_project_summary(project_dir)
  
  # retrieve all variable model parameters
  input_opt_design   <- .rstride$get_variable_model_param(project_summary)
  
  # stop if there are no different input parameter values
  if(nrow(input_opt_design)<=1){
    
    # terminal message
    smd_print('NO VARYING INPUT PARAMETERS IN THE SUMMARY DATA')
    
    # stop
    return(.rstride$no_return_value())
  }
  
  # calculate values and labels for the second y-axis [cases vs. incidence]
  pop_size          <- median(project_summary$num_cases / project_summary$AR)
  
  # OPEN PDF STREAM
  .rstride$create_pdf(project_dir,'summary_inspection',10,7)
  
  # set y-asis ticks and labels
  y_ticks            <- pretty(project_summary$num_cases,10)
  y_ticks_label      <- paste0(round(y_ticks/1e3),'k')
  y_ticks_pop        <- pretty(project_summary$num_cases,5)
  y_ticks_pop_label  <- round(y_ticks_pop/pop_size,digit=2)
  
  # loop over the changing input parameters => plot cases and incidence
  #par(mfrow=c(2,2))
  par(mar = c(10, 4, 4, 4) + 0.3)  # Leave space for 3rd axis
  for(i in 1:ncol(input_opt_design)){
    bxplt <- boxplot(num_cases ~ project_summary[,colnames(input_opt_design)[i]],
            data = project_summary,
            xlab = colnames(input_opt_design)[i],
            ylab = '',cex.axis=0.8,las=2,
            yaxt='n')
    
    # add modified y-axis on the left
    axis(2, at = y_ticks , labels =y_ticks_label,las=2,cex.axis=0.8)
    abline(h=y_ticks,lty=3,col='lightgray')
    mtext("number of cases", side=2, line=2,cex=0.9)
    
    # add y-axis on the right

    axis(4, at = y_ticks_pop , labels =y_ticks_pop_label,cex.axis=0.8)
    mtext("incidence", side=4, line=2,cex=0.9)
    
  }
  
  par(mar = c(10, 4, 1, 4) + 0.3)  # Leave space for 3rd axis
  boxplot(as.formula(paste("num_cases ~ interaction(",paste(colnames(input_opt_design),collapse=','),',drop=T)')),
          data = project_summary,
          ylab = '',
          las=2,
          cex.axis=0.8)
  axis(4, at = y_ticks , labels = y_ticks_label )
  mtext("incidence", side=4, line=2,cex=0.9)
  mtext("total number of cases", side=2, line=2,cex=0.9)
  
  ## SECONDARY CASES ####
  # secondary cases  
  project_summary$num_sec_cases     <- project_summary$num_cases - project_summary$num_infected_seeds
  project_summary$avg_num_sec_cases <- project_summary$num_sec_cases / project_summary$num_infected_seeds
  
  # plot
  sec_cases_formula  <- as.formula(paste('avg_num_sec_cases ~ interaction(',paste(colnames(input_opt_design),collapse=' , '),',drop=TRUE)'))
  bxplot <- boxplot(sec_cases_formula, 
                    data=project_summary,
                    ylab = '',
                    xlab=paste(colnames(input_opt_design)),
                    las=2,
                    cex.axis=0.8)
  mtext("Sec. cases per index case", side=2, line=2,cex=0.9)
  grid()
  
  sec_cases_mean <- aggregate(sec_cases_formula,data=project_summary,mean)
  if(all(sec_cases_mean[,1] == bxplot$names)){
    points(1:nrow(sec_cases_mean),sec_cases_mean$avg_num_sec_cases,pch=4,lwd=2)
    legend('topleft','mean',pch=4,col=1,cex=0.8)
  }
  
  ## SUBSETS ####
  num_subsets <- 10
  project_summary$subset_id <- sample(num_subsets,nrow(project_summary),replace = T)
  sec_cases_formula  <- as.formula(paste('avg_num_sec_cases ~ interaction(subset_id,',paste(colnames(input_opt_design),collapse=' , '),',drop=TRUE)'))
  bxplot <- boxplot(sec_cases_formula, 
                    data=project_summary,
                    ylab = '',
                    xlab=paste(colnames(input_opt_design),collapse=' , '),
                    las=2,
                    cex.axis=0.8)
  mtext("Sec. cases per index case", side=2, line=2,cex=0.9)
  grid()
  
  sec_cases_mean <- aggregate(sec_cases_formula,data=project_summary,mean)
  points(which(!is.na(bxplot$stats[1,])),sec_cases_mean$avg_num_sec_cases,pch=4,lwd=2)
  abline(v=(0:20)*num_subsets+0.5)
  legend('topleft','mean',pch=4,col=1,cex=0.8)
  title(paste('subset size:',mean(bxplot$n)))
  
  dev.off()
  
  # CLI message
  smd_print('INSPECTION OF SUMMARY DATA COMPLETE')
  
}

## HELP FUNCTION ####
.rstride$get_variable_model_param <- function(project_summary){
  
  if(ncol(project_summary) == 1){
    return(unique(project_summary))
  }
  
  input_opt    <- .rstride$get_unique_param_list(project_summary)
  input_opt    <- input_opt[lapply(input_opt,length)>1]
  
  # get parameter combinations
  input_opt_design <- unique(project_summary[,names(input_opt)])
  
  # with only one parameter, convert vector into nx2 matrix with dummy column
  if(length(input_opt)==1){
    input_opt_design <- data.frame(input_opt,
                                   num_days = unique(project_summary$num_days) )
  }
  
  # with only identical parameters, use the r0
  if(length(input_opt)==0){
    input_opt_design <- data.frame(r0=unique(project_summary$r0),
                                   num_days=unique(project_summary$num_days))
  }
  
  return(input_opt_design)
  
}

## HELP FUNCTION ####
.rstride$get_unique_param_list <- function(project_summary){

  # option to select "holidays_calendar" basic file name by excluding the relative path, to overcome (experiment-specific) file names 
  if('holidays_file' %in% names(project_summary)){
    project_summary$holidays_file <- basename(project_summary$holidays_file)
  }
  
  col_output <- c('run_time', 'total_time', 'num_cases', 'AR' )
  col_extra  <- c('rng_seed','output_prefix','transmission_probability','exp_id','config_id','contact_id')
  col_poison <- names(project_summary)[grepl('_pois',names(project_summary))]
  col_input  <- !(names(project_summary) %in% c(col_output,col_extra,col_poison))
  
  # get unique values per parameter
  input_opt    <- lapply(project_summary[,col_input],unique)
  
  # return
  return(input_opt)
}

## HELP FUNCTION ####
.rstride$get_config_id <- function(project_matrix){
    
  # retrieve all variable model parameters
  input_opt_design   <- .rstride$get_variable_model_param(project_matrix)
  
  # collapse variables of input_opt_design (or just copy if only one column)
  if(ncol(input_opt_design) == 1) {
    config_id  <- project_matrix[,colnames(input_opt_design)]
  } else{
    config_id  <- apply(project_matrix[,names(input_opt_design)],1,paste, collapse='_')
  }
  
  return(config_id)
}

.rstride$get_contact_id <- function(project_matrix){
  
  # retrieve all variable model parameters
  input_opt_design   <- .rstride$get_variable_model_param(project_matrix)
  
  # add contact id
  flag_opt_input_tracing <- (grepl('cnt_reduction',names(input_opt_design)) & grepl('exit',names(input_opt_design)))
  if(any(flag_opt_input_tracing)){
    colnames_contact           <- names(input_opt_design)[flag_opt_input_tracing]
    if(length(colnames_contact) == 1) {
      contact_id  <- project_matrix[,colnames_contact]
    } else{
      contact_id  <- apply((1-project_matrix[,colnames_contact])*100,1,paste, collapse='_')
    }
  } else{
    contact_id  <- .rstride$get_config_id(project_matrix)
  }

  return(contact_id)
}

