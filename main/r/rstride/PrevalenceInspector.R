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
  data_prevalence               <- .rstride$load_aggregated_output(project_dir,'data_prevalence')
  
  if(length(data_prevalence) == 1 && is.na(data_prevalence)){
    smd_print('NO PREVALENCE DATA AVAILABLE.')
    return(NA)
  }
  
  # get specific prevalence output
  data_prevalence_infected    <- get_prevalence_matrix(data_prevalence,'prevalence_infected')
  data_prevalence_exposed     <- get_prevalence_matrix(data_prevalence,'prevalence_exposed')
  data_prevalence_infectious  <- get_prevalence_matrix(data_prevalence,'prevalence_infectious')
  data_prevalence_symptomatic <- get_prevalence_matrix(data_prevalence,'prevalence_symptomatic')
  
  #TODO: re-use Incidence 'pcolor'
  # set color definitions and other layout definitions
  pcolor <- data.frame(E = "black",  # exposed (or total infections)
                       I = "darkgoldenrod3",  # infectious
                       S = "red",  # symptomatic
                       alpha = 0.1,
                       lwd = 3,
                       stringsAsFactors = F)  # data
  
  # change transparancey for low number of rng-runs
  if(max(table(project_summary$exp_id)) < 10){
    pcolor$alpha <- 0.2
  }
  
  .rstride$create_pdf(project_dir,'prevalence',width = 6, height = 2.5)
  par(mar=c(3,5,1,3))
  
  # col_days <- which(grepl('day',names(data_prevalence_exposed)))
  sim_dates_all <- range(as.Date(project_summary$start_date[1]) + (1:ncol(data_prevalence))-1)
  #sim_dates[2] <- sim_dates[2] - 68
  
  sim_dates_opt <- list(sim_dates_all,
                        c(as.Date("2020-03-01"),as.Date("2020-03-15")),
                        c(as.Date("2020-05-18"),as.Date("2020-06-15")),
                        c(as.Date(project_summary$start_date[1]),as.Date("2020-06-01")))
  
  for(sim_dates in sim_dates_opt){
    
    y_lim <- range(data_prevalence_exposed,data_prevalence_infectious,na.rm = T)
    if(!any(sim_dates == sim_dates_all[1])){
      y_lim <- c(0,3e4)
    }

  plot(sim_dates,y_lim,
       col=0,
       xlab='Time',
       ylab='Prevalence',
       xaxt='n',
       yaxt='n')
  add_x_axis(sim_dates)
  add_y_axis(y_lim)
  add_breakpoints()
  add_legend_prevalence(pcolor,'topright')
  
  i_exp <- 1
  for(i_exp in 1:nrow(data_prevalence_exposed)){
    
    flag_exp <- project_summary$exp_id == i_exp
    sim_date <- as.Date(project_summary$start_date[flag_exp]) + (1:ncol(data_prevalence_exposed))-1
      lines(x = sim_date,
           y = data_prevalence_exposed[i_exp,],
           col = pcolor$E)
      lines(x = sim_date,
            y = data_prevalence_infectious[i_exp,],
            col = pcolor$I)
      lines(x = sim_date,
            y = data_prevalence_symptomatic[i_exp,],
            col = pcolor$S)
  }
  
  }
  
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

# define the legend with all categories
add_legend_prevalence <- function(pcolor,legend_pos = 'topleft'){
  
  legend(legend_pos,
         c('Exposed (latent)',
           'Infectious',
           'Symptomatic'),
         col=unlist(pcolor),
         lwd=c(2,2,2),
         cex=0.5,
         bg='white')
}


