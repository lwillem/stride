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
#  Copyright 2020, Willem L, Kuylen E & Broeckhove J
############################################################################# #
#
# MODEL CONTACT TRACING EXPLORATION
#
############################################################################# #

#' @param project_dir   name of the project folder
#' @param num_selection the number of experiments with minimal LS score to select and present
inspect_tracing_data <- function(project_dir)
{
  # command line message
  smd_print('INSPECT CONTACT TRACING DATA...')
  
  # load project summary
  project_summary    <- .rstride$load_project_summary(project_dir)
  
  # retrieve all variable model parameters
  input_opt_design   <- .rstride$get_variable_model_param(project_summary)
  
  # get all tracing output
  data_tracing_all <- .rstride$load_aggregated_output(project_dir,'data_tracing')
  
  if(length(data_tracing_all) == 1 && is.na(data_tracing_all)){
    smd_print('NO CONTACT TRACE DATA AVAILABLE.')
    return(NA)
  }
  
 
  # add config_id 
  # get variable names of input_opt_design (fix if only one column)
  if(ncol(input_opt_design) == 1) {
    project_summary$config_id  <- project_summary[,colnames(input_opt_design)]
    input_opt_design           <- data.frame(input_opt_design,config_id = c(input_opt_design))
  } else{
    project_summary$config_id  <- apply(project_summary[,names(input_opt_design)],1,paste, collapse='_')
    input_opt_design$config_id <- apply(input_opt_design,1,paste, collapse='_')
  }
  
  # add config_id to incidence data
  data_tracing_all         <- merge(data_tracing_all,project_summary[,c('exp_id','config_id','start_date')] )

  ## ENSEMBLE  ####
  .rstride$create_pdf(project_dir,'contact_tracing_all',width = 4, height = 4)
  par(mar=c(6,5,4,1))
  
  opt_config <- unique(data_tracing_all$config_id)
  i_config <- 2
  for(i_config in 1:length(opt_config)){
    print(i_config)
    data_tracing_sel <- data_tracing_all[data_tracing_all$config_id == opt_config[i_config],]
    
    # add date
    sim_start_date <- as.Date(unique(data_tracing_sel$start_date))
    data_tracing_sel$sim_day_date <- sim_start_date + data_tracing_sel$sim_day
    
    # add value to aggregate
    data_tracing_sel$num_tests <- 1
    head(data_tracing_sel)
    
    # index cases
    data_tracing_index       <- data_tracing_sel[data_tracing_sel$pool_type == '-1',]
    tracing_num_day_index    <- aggregate(num_tests ~ sim_day_date + sim_day + exp_id + config_id, data = data_tracing_index, sum)
    tracing_num_day_contacts <- aggregate(num_contacts_tested ~ sim_day_date + exp_id + config_id, data = data_tracing_index, sum)
    tracing_num_day_contacts_mean <- aggregate(num_contacts_tested ~ sim_day_date + exp_id + config_id, data = data_tracing_index, mean)
    
    
    # contacts in quarentine
    data_tracing_contacts      <- data_tracing_sel[data_tracing_sel$pool_type != '-1',]
    tracing_num_day_identified <- aggregate(num_tests ~ sim_day_date + exp_id + config_id, data = data_tracing_contacts, sum)
    tracing_num_day_sympt      <- aggregate(num_tests ~ sim_day_date + sim_day + exp_id + config_id + is_symptomatic, data = data_tracing_contacts, sum)
    tracing_num_day_sympt_mean <- aggregate(num_tests ~ sim_day_date + config_id + is_symptomatic, data = tracing_num_day_sympt, mean)
    
    
    y_lim <- range(0,tracing_num_day_index$num_tests)
    boxplot(num_tests ~ sim_day_date, data=tracing_num_day_index,
            ylab='index cases',main=opt_config[i_config],
            ylim = y_lim,las=2,xlab='',cex=0.8,
            xaxt='n')
    
    x_ticks_label <- format(unique(tracing_num_day_index$sim_day_date),'%d/%m')
    x_ticks       <- seq(1,length(x_ticks_label),7)
    x_ticks_label <- x_ticks_label[x_ticks]
    axis(1,x_ticks,x_ticks_label,las=2)
    grid(nx=NA,ny=NULL)
    abline(v=x_ticks,lty=3,col='lightgray')
    
    
    plot(tracing_num_day_sympt_mean$sim_day_date,
         tracing_num_day_sympt_mean$num_tests,
         col=tracing_num_day_sympt_mean$is_symptomatic+1,
         type='p',
         pch=16,
         xlab='',
         ylab='Secondary cases identified and isolated',
         las=2,
         main=opt_config[i_config])
    legend('topleft',
           c('asymptomatic',
             'symptomatic'),
           pch=16,
           col=1:2,
           cex=0.7)
    
    boxplot(num_contacts_tested ~ sim_day_date, data=tracing_num_day_contacts,
            las=2,ylab='total number of contacts tested',main=opt_config[i_config],
            xaxt='n')
    
    x_ticks_label <- format(unique(tracing_num_day_contacts$sim_day_date),'%d/%m')
    x_ticks       <- seq(1,length(x_ticks_label),7)
    x_ticks_label <- x_ticks_label[x_ticks]
    axis(1,x_ticks,x_ticks_label,las=2)
    grid(nx=NA,ny=NULL)
    abline(v=x_ticks,lty=3,col='lightgray')
    
    
    boxplot(num_contacts_tested ~ sim_day_date, data = data_tracing_index,
            las=2,ylab='contacts tested per index case',main=opt_config[i_config],
            xaxt='n')
    
    x_ticks_label <- format(unique(data_tracing_index$sim_day_date),'%d/%m')
    x_ticks       <- seq(1,length(x_ticks_label),7)
    x_ticks_label <- x_ticks_label[x_ticks]
    axis(1,x_ticks,x_ticks_label,las=2)
    grid(nx=NA,ny=NULL)
    abline(v=x_ticks,lty=3,col='lightgray')

    boxplot(num_contacts_tested ~ sim_day_date, data = data_tracing_index,outline=F,
            las=2,ylab='contacts tested per index case',main=opt_config[i_config],
            xaxt='n')
    
    axis(1,x_ticks,x_ticks_label,las=2)
    grid(nx=NA,ny=NULL)
    abline(v=x_ticks,lty=3,col='lightgray')
    
    
    }
  
  # close pdf stream
  dev.off()
  
    # command line message
  smd_print('INSPECTION OF CONTACT TRACING DATA COMPLETE')
}
