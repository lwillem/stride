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
    smd_print('NO CONTACT TRACING DATA AVAILABLE.')
    return(.rstride$no_return_value())
  }
  
  # make sure that almost all values are stored as integers
  pool_type_col <- names(data_tracing_all) %in% c('pool_type')
  data_tracing_all[,!pool_type_col] <- data.frame(apply(data_tracing_all[,!pool_type_col], 2, as.integer))
  dim(data_tracing_all)

  # add config_id 
  project_summary$config_id <- .rstride$get_config_id(project_summary)
  
  # add config_id to incidence data
  data_tracing_all         <- merge(data_tracing_all,project_summary[,c('exp_id','config_id','start_date')] )

  dim(data_tracing_all)
  length(unique(data_tracing_all$case_id))
  print(table(data_tracing_all$pool_type))
  
  ## ENSEMBLE  ####
  .rstride$create_pdf(project_dir,'contact_tracing_all',width = 14, height = 8)
  par(mar=c(6,5,4,1))
  
  opt_config <- unique(data_tracing_all$config_id)
  i_config <- 1
  for(i_config in 1:length(opt_config)){

    # get subset
    data_tracing_sel <- data_tracing_all[data_tracing_all$config_id == opt_config[i_config],]
    
    # add date
    sim_start_date <- as.Date(unique(data_tracing_sel$start_date))
    data_tracing_sel$sim_day_date <- sim_start_date + data_tracing_sel$sim_day
    
    # add value to aggregate
    data_tracing_sel$num_tests <- 1
    head(data_tracing_sel)
    
    # index cases
    data_tracing_index       <- data_tracing_sel[data_tracing_sel$pool_type == 'Index',]
    tracing_num_day_index    <- aggregate(num_tests ~ sim_day_date + sim_day + exp_id + config_id, data = data_tracing_index, sum)
    tracing_num_day_contacts <- aggregate(num_contacts_tested ~ sim_day_date + exp_id + config_id, data = data_tracing_index, sum)
    tracing_num_day_contacts_mean <- aggregate(num_contacts_tested ~ sim_day_date + exp_id + config_id, data = data_tracing_index, mean)
    
    
    # contacts in quarentine
    data_tracing_contacts      <- data_tracing_sel[data_tracing_sel$pool_type != 'Index',]
    tracing_num_day_identified <- aggregate(num_tests ~ sim_day_date + exp_id + config_id, data = data_tracing_contacts, sum)
    tracing_num_day_sympt      <- aggregate(num_tests ~ sim_day_date + sim_day + exp_id + config_id + is_symptomatic, data = data_tracing_contacts, sum)
    tracing_num_day_sympt_mean <- aggregate(num_tests ~ sim_day_date + config_id + is_symptomatic, data = tracing_num_day_sympt, mean)
    
    
    y_lim <- range(0,tracing_num_day_index$num_tests)
    boxplot(num_tests ~ sim_day_date, data=tracing_num_day_index,
            ylab='index cases',main=opt_config[i_config],
            ylim = y_lim,las=2,xlab='',cex=0.8,
            xaxt='n',yaxt='n')
    add_y_axis(tracing_num_day_index$num_tests)

    x_ticks_label <- format(unique(tracing_num_day_index$sim_day_date),'%d/%m')
    x_ticks       <- seq(1,length(x_ticks_label),14)
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
         main=opt_config[i_config],
         xaxt='n')
    add_x_axis(tracing_num_day_sympt_mean$sim_day_date)
    legend('topleft',
           c('asymptomatic',
             'symptomatic'),
           pch=16,
           col=1:2,
           cex=0.7)
    
    boxplot(num_contacts_tested ~ sim_day_date, data=tracing_num_day_contacts,
            las=2,ylab='total number of contacts tested',main=opt_config[i_config],
            xaxt='n',yaxt='n')
    add_y_axis(tracing_num_day_contacts$num_contacts_tested)
    x_ticks_label <- format(unique(tracing_num_day_contacts$sim_day_date),'%d/%m')
    x_ticks       <- seq(1,length(x_ticks_label),14)
    x_ticks_label <- x_ticks_label[x_ticks]
    
    axis(1,x_ticks,x_ticks_label,las=2)
    grid(nx=NA,ny=NULL)
    abline(v=x_ticks,lty=3,col='lightgray')
    
    
    boxplot(num_contacts_tested ~ sim_day_date, data = data_tracing_index,
            las=2,ylab='contacts tested per index case',main=opt_config[i_config],
            xaxt='n',yaxt='n')
    add_y_axis(data_tracing_index$num_contacts_tested)
    x_ticks_label <- format(unique(data_tracing_index$sim_day_date),'%d/%m')
    x_ticks       <- seq(1,length(x_ticks_label),7)
    x_ticks_label <- x_ticks_label[x_ticks]
    axis(1,x_ticks,x_ticks_label,las=2)
    grid(nx=NA,ny=NULL)
    abline(v=x_ticks,lty=3,col='lightgray')
    
    boxplot(num_unique_contacts ~ sim_day_date, data = data_tracing_index,
            las=2,ylab='unique contacts per index case',main=opt_config[i_config],
            xaxt='n',yaxt='n')
    add_y_axis(data_tracing_index$num_contacts_tested)
    x_ticks_label <- format(unique(data_tracing_index$sim_day_date),'%d/%m')
    x_ticks       <- seq(1,length(x_ticks_label),7)
    x_ticks_label <- x_ticks_label[x_ticks]
    axis(1,x_ticks,x_ticks_label,las=2)
    grid(nx=NA,ny=NULL)
    abline(v=x_ticks,lty=3,col='lightgray')

    
    ## LOCATION
    data_tracing_contacts$cases <- 1
    dt_location <- aggregate(cases ~ pool_type + exp_id + sim_day_date, data = data_tracing_contacts, sum )
    dt_location_mean <- aggregate(cases ~ pool_type + sim_day_date, data = dt_location, mean )
   
    table(data_tracing_contacts$pool_type) / nrow(data_tracing_contacts)
    
    plot(cases ~ sim_day_date, 
         data = dt_location_mean[dt_location_mean$pool_type == 'Household',],
         ylim = range(dt_location_mean$cases),xaxt='n',yaxt='n',
         xlab='',ylab='Secondary cases identified and isolated')
    add_x_axis(dt_location_mean$sim_day_date)
    add_y_axis(dt_location_mean$cases)
    
    points(cases ~ sim_day_date, 
         data = dt_location_mean[dt_location_mean$pool_type == 'Workplace',],
         col = 2)
    points(cases ~ sim_day_date, 
         data = dt_location_mean[dt_location_mean$pool_type == 'Community',],
         col = 3)
    points(cases ~ sim_day_date, 
           data = dt_location_mean[dt_location_mean$pool_type == 'School',],
           col = 4)
    legend('topleft',
           c('Household','Workplace','Community','School'),
           fill = 1:4,
           cex=0.5)
    }
  
  # close pdf stream
  dev.off()
  
    # command line message
  smd_print('INSPECTION OF CONTACT TRACING DATA COMPLETE')
}

