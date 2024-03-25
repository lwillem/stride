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
# MODEL PREVALENCE EXPLORATION
#
############################################################################# #

#' @param project_dir   name of the project folder
inspect_prevalence_data <- function(project_dir)
{
  # command line message
  smd_print('INSPECT PREVALENCE DATA...')
  
  # load project summary
  project_summary    <- .rstride$load_project_summary(project_dir)
  
  # retrieve all variable model parameters
  input_opt_design   <- .rstride$get_variable_model_param(project_summary)
  
  # get all prevalence output
  data_prevalence_all               <- .rstride$load_aggregated_output(project_dir,'data_prevalence')
  
  if(length(data_prevalence_all) == 1 && is.na(data_prevalence_all)){
    smd_print('NO PREVALENCE DATA AVAILABLE.')
    return(NA)
  }
  
  #TODO: re-use Incidence 'pcolor'
  # set color definitions and other layout definitions
  pcolor <- data.frame(E = "black",  # exposed (or total infections)
                       I = "darkgoldenrod3",  # infectious
                       S = "red",  # symptomatic
                       H = "blue",
                       alpha = 0.1,
                       lwd = 3,
                       stringsAsFactors = F)  # data
  
  # change transparancey for low number of rng-runs
  if(max(table(project_summary$exp_id)) < 10){
    pcolor$alpha <- 0.2
  }
  
  # open pdf stream
  .rstride$create_pdf(project_dir,'prevalence',width = 14, height = 8)
  par(mar=c(3,5,1,3))
  
  i_config <- 1
  for(i_config in 1:nrow(input_opt_design)){
  
    # subset transmission output corresponding the 'input_opt_design' row
    flag_exp            <- .rstride$get_equal_rows(project_summary,input_opt_design[i_config,])
    data_prevalence     <- data_prevalence_all[data_prevalence_all$exp_id %in% project_summary$exp_id[flag_exp],]
    dim(data_prevalence)
  
    if(nrow(data_prevalence) > 0)
    {
      # get specific prevalence output
      data_prevalence_infected     <- get_prevalence_matrix(data_prevalence,'prevalence_infected')
      data_prevalence_exposed      <- get_prevalence_matrix(data_prevalence,'prevalence_exposed')
      data_prevalence_infectious   <- get_prevalence_matrix(data_prevalence,'prevalence_infectious')
      data_prevalence_symptomatic  <- get_prevalence_matrix(data_prevalence,'prevalence_symptomatic')
      data_prevalence_hospitalized <- get_prevalence_matrix(data_prevalence,'prevalence_hospitalised')
      data_prevalence_date         <- get_prevalence_dates(data_prevalence)
      
      sim_dates <- range(data_prevalence_date)
      y_lim     <- range(0,data_prevalence_exposed,data_prevalence_infectious,na.rm = T)
      
      plot(sim_dates,y_lim,
           col=0,
           xlab='Time',
           ylab='Prevalence',
           xaxt='n',
           yaxt='n')
      add_x_axis(sim_dates)
      add_y_axis(y_lim)
      add_intervention_dates(project_summary[flag_exp,])
      add_legend_prevalence(pcolor,'topright')
      
      i_exp <- 1
      for(i_exp in 1:nrow(data_prevalence_exposed)){
        
          lines(x = data_prevalence_date,
               y = data_prevalence_exposed[i_exp,],
               col = pcolor$E)
          lines(x = data_prevalence_date,
                y = data_prevalence_infectious[i_exp,],
                col = pcolor$I)
          lines(x = data_prevalence_date,
                y = data_prevalence_symptomatic[i_exp,],
                col = pcolor$S)
          lines(x = data_prevalence_date,
                y = data_prevalence_hospitalized[i_exp,],
                col = pcolor$H)
      } # end if-clause, nrow(data_prevalence)>0 
    } # end for-loop i_exp
  } # end for-loop config_id
  
  # close pdf stream
  dev.off()
}

get_prevalence_matrix <- function(data_prevalence,col_name){
  num_days <- length(unique(data_prevalence$sim_date))
  num_exp  <- length(unique(data_prevalence$exp_id))
  return(matrix(data_prevalence[,col_name],
                nrow=num_exp,
                ncol=num_days,
                byrow = TRUE))
}

get_prevalence_dates <- function(data_prevalence){
  return(as.Date(unique(data_prevalence$sim_date)))
}

# define the legend with all categories
add_legend_prevalence <- function(pcolor,legend_pos = 'topleft'){
  
  legend(legend_pos,
         c('Exposed (latent)',
           'Infectious',
           'Symptomatic',
           'Hospitalized'),
         col=unlist(pcolor),
         lwd=2,
         cex=0.5,
         bg='white')
}


