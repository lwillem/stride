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
#  Copyright 2020, Willem L
############################################################################ #

############################################################################# #
# EXPLORE PARTICIPANT DATA                                                 ####
############################################################################# #

inspect_participant_data <- function(project_dir, save_pdf = TRUE)
{
  # command line message
  smd_print('INSPECT PARTICIPANT DATA...')
  
  # load project summary
  project_summary <- .rstride$load_project_summary(project_dir)
  
  # if no participants are surveyed, stop
  if(any(project_summary$num_participants_survey==0)){
    # command line message
    smd_print('NO PARTICIPANTS SELECTED')
    return(.rstride$no_return_value())
  }
  
  # get all transmission output
  data_participants_all      <- .rstride$load_aggregated_output(project_dir,'data_participants')
  if(all(is.na(data_participants_all))){
    # command line message
    smd_print('NO PARTICIPANT DATA AVAILABLE')
    return(.rstride$no_return_value())
  }
  
  # retrieve all variable model parameters
  input_opt_design        <- .rstride$get_variable_model_param(project_summary)
  
  # open pdf stream
  if(save_pdf) .rstride$create_pdf(project_dir,'survey_participant_inspection',10,7)

  # set maximum number of days symtomatic/infectious 
  max_days_infection <- max(c(30,unlist(data_participants_all[,c("end_infectiousness","end_symptomatic")])),na.rm=T)
  
  i_config <- 1
  for(i_config in 1:nrow(input_opt_design))
  {
    # (re)set figure panels
    if(save_pdf) par(mfrow=c(2,4))
    
    # select the participant output subset, corresponding the 'input_opt_design' row
    flag_exp            <- .rstride$get_equal_rows(project_summary,input_opt_design[i_config,])
    data_part           <- .rstride$load_aggregated_output(project_dir,'data_participants',project_summary$exp_id[flag_exp])
    num_runs_exp        <- sum(flag_exp)

    # check if participant data is avaiable, if not, go to next iteration
    if(nrow(data_part)==0){
      next
    }
    
    # adjust type of column
    data_part$start_symptomatic  <- as.numeric(data_part$start_symptomatic)
    data_part$end_infectiousness <- as.numeric(data_part$end_infectiousness)
    
    # adjust for asymptomatic cases
    flag_asymptomatic <- data_part$start_symptomatic == data_part$end_symptomatic
    data_part$start_symptomatic[flag_asymptomatic] <- NA
    data_part$end_symptomatic[flag_asymptomatic] <- NA
    
    num_part <- nrow(data_part)
    freq_start_inf  <- table(data_part$start_infectiousness) / num_part
    freq_start_symp <- table(data_part$start_symptomatic)    / num_part
    data_part[1,]
    all_inf  <- matrix(0,num_part,max_days_infection)
    all_symp <- matrix(0,num_part,max_days_infection)
    #note: if start == end, no infectious/symptomatic stage has been present
    for(i in 1:num_part){
      all_inf[i,data_part$start_infectiousness[i]:data_part$end_infectiousness[i]] <- 1
      if(!is.na(data_part$start_symptomatic[i]))
      all_symp[i,data_part$start_symptomatic[i]:data_part$end_symptomatic[i]] <- 1
    }

    plot(1:max_days_infection,colMeans(all_inf),
         ylab='population fraction',
         xlab='days since infection',
         type='b',lwd=3,col=2,
         ylim=c(0,1))
    points(1:max_days_infection,colMeans(all_symp),lwd=3,col=4,type='b')
    legend('topright',c('infectious','symptomatic'),col=c(2,4),lwd=4,cex=0.8)
    abline(v=6:9,lty=3)
    
    f_data <- data_part$start_symptomatic; f_main <- 'debug'
    plot_cum_distr <- function(f_data,f_main,f_x_lab = 'period (days)'){
      tbl_data <- table(f_data)/length(f_data)
      tbl_data_cumm <- cumsum(tbl_data)
      plot(tbl_data,xlim=c(0,20),ylim=0:1,
           xlab=f_x_lab,ylab='frequency',main=f_main)
      lines(as.numeric(names(tbl_data_cumm)),tbl_data_cumm,col=4,lwd=2,type='b')
      
      legend_position <- 'topleft'
      if(max(as.numeric(names(tbl_data)))<10) {legend_position <- 'topright'}
      legend(legend_position,c('per day','cumulative'),col=c(1,4),lwd=2,cex=0.8)
      grid()
      abline(h=0.5)
    }
    
    # days asymptomitic & infectious ####
    start_symtomatic <- data_part$start_symptomatic
    start_symtomatic[is.na(data_part$start_symptomatic)] <- data_part$end_infectiousness[is.na(data_part$start_symptomatic)]
    days_asymptomatic_infectious <- start_symtomatic - data_part$start_infectiousness
   
    ## distributions ####   
    plot_cum_distr(data_part$start_infectiousness,f_main='start_infectiousness',f_x_lab='time since infection (days)')
    plot_cum_distr(data_part$start_symptomatic,f_main='start_symptomatic',f_x_lab='time since infection (days)')
    plot_cum_distr(days_asymptomatic_infectious,f_main='days infectious \n& not symptomatic')
    plot_cum_distr(data_part$end_infectiousness-data_part$start_infectiousness,f_main='days infectious')
    plot_cum_distr(data_part$end_symptomatic-data_part$start_symptomatic,f_main='days symptomatic')
    
    # fraction symptomatic by age ####
    tbl_sympt_age <- table(!is.na(data_part$start_symptomatic),data_part$part_age)
    plot(tbl_sympt_age[2,]/colSums(tbl_sympt_age),
         xlab='age',
         ylab = 'Relative frequency symptomatic')
    abline(h=sum(tbl_sympt_age[2,]) / sum(tbl_sympt_age))
    legend('topleft',
           'mean',
           col = 1,
           lwd = 1,
           cex=0.8)
    
    # ## POPULATION
    population_age <- as.data.frame(table(part_age = data_part$part_age))
    # 
    # ## POPULATION IMMUNITY
    # data_part$is_immune   <- data_part$is_immune == "TRUE"
    # 
    # immune_age <- data.frame(table(is_immune = data_part$is_immune, part_age = data_part$part_age),stringsAsFactors = F)
    # immune_age$part_age <- as.numeric(levels(immune_age$part_age)[(immune_age$part_age)])
    # names(immune_age)  
    # flag <- immune_age$is_immune == FALSE
    # plot(immune_age$part_age[flag],
    #      immune_age$Freq[flag]/population_age$Freq,
    #      xlab='age',
    #      ylab='fraction susceptible',
    #      main='population susceptibility',
    #      pch=19, lwd=3, ylim=0:1
    # )
    # 
    # names(data_part)
    
   
    ## SCHOOLING ----
    data_part$enrolled_school <- data_part$school_id != 0 | data_part$college_id != 0
    school_age <- data.frame(table(school_enrolled = data_part$enrolled_school,part_age = data_part$part_age))
    school_age$part_age <- as.numeric(levels(school_age$part_age)[(school_age$part_age)])
    flag <- school_age$school_enrolled == TRUE
    plot(school_age$part_age[flag],
         school_age$Freq[flag]/population_age$Freq,
         xlab='age',
         ylab='population fraction',
         main='population enrolled in school',
         pch=1, lwd=3, xlim=c(0,30)
    )
    abline(v=c(0,3,6,12,18,25)-0.5)
    text(x=0,y=0.02,'kindergarten',srt=90,pos=4)
    text(x=3,y=0.02,'pre-school',srt=90,pos=4)
    text(x=7,y=0.02,'primary school',srt=90,pos=4)
    text(x=13,y=0.02,'secundary school',srt=90,pos=4)
    text(x=20,y=0.02,'tertiary school',srt=90,pos=4)
    # 
    
    # tmp_age      <- table(data_part$part_age)
    # tmp_college  <- table(data_part$part_age,data_part$college_id != 0)
    # tmp_school   <- table(data_part$part_age,data_part$school_id != 0)
    # 
    # barplot(rbind(tmp_school[,2]/tmp_age,tmp_college[,2]/tmp_age),
    #         xlab='age',
    #         ylab='relative fraction')
    # legend('topright',
    #        c('school open',
    #          'school-closed'),
    #        fill=grey.colors(2))   
    
    
  }
  
  # close PDF stream
  if(save_pdf) dev.off()
  
  # command line message
  smd_print('INSPECTION OF PARTICIPANT DATA COMPLETE')
  
} # function end
