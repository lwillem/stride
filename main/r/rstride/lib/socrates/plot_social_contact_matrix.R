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
#  Copyright 2024
############################################################################ #
#
# PLOT SOCIAL CONTACT DATA
#___________________________________________________________________________

# Loading packages
socrates_packages <- c('socialmixr','countrycode','data.table')
smd_load_packages(socrates_packages)

# requires 'simage' and 'splot' functions from npsp package (removed from CRAN)

# plot matrices according the Socrates app
plot_socrates_all <- function(data_cnt,data_part,age_cat_breaks,project_dir,exp_tag,survey_start){
  
  .rstride$create_pdf(project_dir,paste0(exp_tag,'_cnt_matrix_all'),9,4)
  
  opt_day <- unique(data_cnt$sim_day)
  i_day <- 1
  for(i_day in opt_day){
    data_cnt_day <- data_cnt[data_cnt$sim_day == i_day,]
    plot_socrates_location(data_cnt_day,data_part,age_cat_breaks,as.Date(survey_start) + i_day)
  }
  
  dev.off()
  
  if(any(data_part$start_symptomatic %in% opt_day))
  {
    .rstride$create_pdf(project_dir,paste0(exp_tag,'_cnt_matrix_symptomatic'))
    # symptomatic
    for(i_day in opt_day){
     
      # select contacts of 'i_day' and symptomatic participants (of the infected seeds)
      data_part_sympt <- data_part[data_part$is_infected == 1  & (data_part$start_symptomatic-1) <= i_day & (data_part$end_symptomatic-1) > i_day,]
      data_cnt_day    <- data_cnt[data_cnt$sim_day == i_day & data_cnt$local_id %in% data_part_sympt$local_id,]
      
      if(nrow(data_part_sympt)>0 & nrow(data_cnt_day)>0){
        print(i_day)  
        plot_socrates_location(data_cnt_day,data_part_sympt,age_cat_breaks,as.Date(survey_start) + i_day)
      }
    }
     dev.off()
  } # end if-clause symptomatic infections?
}

 # data_cnt <- data_cnt_day; data_part <- data_part_sympt;survey_day <- as.Date(survey_start) + i_day
plot_socrates_location <- function(data_cnt,data_part,age_cat_breaks,survey_day){
  
  par(mfrow=c(2,3))
  
  ## TOTAL
  plot_contact_matrix_socrates(data_cnt,data_part,'total',age_cat_breaks)
  
  ## HOUSEHOLD
  plot_contact_matrix_socrates(data_cnt[data_cnt$cnt_home==1,],data_part,'@household',age_cat_breaks)
  
  ## SCHOOL
  plot_contact_matrix_socrates(data_cnt[data_cnt$cnt_school==1,],data_part[data_part$student==T,],'@school',age_cat_breaks)
  
  ## WORK
  plot_contact_matrix_socrates(data_cnt[data_cnt$cnt_work==1,],data_part[data_part$employed==T,],'@work',age_cat_breaks)
  
  ## PRIMARY COMMUNITY
  plot_contact_matrix_socrates(data_cnt[data_cnt$cnt_prim_comm==1,],data_part,'@weekend community',age_cat_breaks)
  
  ## SECONDARY COMMUNITY
  plot_contact_matrix_socrates(data_cnt[data_cnt$cnt_sec_comm==1,],data_part,'@week community',age_cat_breaks)
  
  ## HOUSEHOLD CLUSTER
  plot_contact_matrix_socrates(data_cnt[data_cnt$cnt_hh_cluster==1,],data_part,'@household cluster',age_cat_breaks)
  
  plot(0,0,col=0,axes=F,xlab='',ylab='')
  # survey_day <- as.Date(exp_summary$start_date) + unique(data_cnt$sim_day)
  survey_day <- as.Date(survey_day)
  text(0,0,paste(format(survey_day,'%A'),survey_day),pos=3)
  text(0,0,paste('Contacts when symptomatic:',sum(data_cnt$part_sympt),
                  '\nNumber of participants:',nrow(data_part)),pos=1)
  
}

# data_cnt <- data_cnt[data_cnt$cnt_work==1,]; data_part <- data_part[data_part$employed==T,]
plot_contact_matrix_socrates <- function(data_cnt,data_part,figure_title,age_cat_breaks){
  
  
  if(nrow(data_cnt)>0 & nrow(data_part)>0){
    
    # get socialmixr 'participants' object
    db_participants   <- data.frame(part_id     = data_part$local_id,
                                    part_age    = data_part$part_age,
                                    part_gender = NA,
                                    country     = "Belgium",
                                    day         = NA,
                                    month       = NA,
                                    year        = 2020,
                                    dayofweek   = NA,
                                    holiday     = FALSE,
                                    weekday     = NA,
                                    stringsAsFactors = F)
    
    # get socialmixr 'contacts' object
    db_contacts       <- data.frame(part_id         = data_cnt$local_id,
                                    cnt_age_exact   = as.integer(round(data_cnt$cnt_age)),
                                    cnt_age_est_min = as.integer(round(data_cnt$cnt_age)),
                                    cnt_age_est_max = as.integer(round(data_cnt$cnt_age)),
                                    data_cnt[,c("cnt_home","cnt_work","cnt_school",
                                                "cnt_prim_comm","cnt_sec_comm","part_sympt","cnt_sympt" )])
    
    # get socialmixr 'survey' object
    survey_rstride <- survey(participants = db_participants,
                             contacts     = db_contacts)
  
    # get matrix
    cnt_matrix <- contact_matrix(survey_rstride,age.limits = age_cat_breaks)  
    
    # plot matrix
    plot_cnt_matrix(cnt_matrix$matrix,figure_title)
  }
}

#mij <- contact_matrix(polymod, countries = "United Kingdom", age.limits = c(0, 1, 5, 15))$matrix
#mij <- matrix_out$matrix
plot_cnt_matrix <- function(mij,plot_title_extra = ''){
  
  if(all(is.na(mij))){
    return(NA)
  }
  
  # set digits
  format_num_digits <- 2
  
  redc <- rev(heat.colors(100))
  par(mar=c(5, 6, 2, 2),mgp=c(3,0.5,0))
  p <- simage(s = mij, 
             xlab="Age of participant (year)",
             ylab="Age of contact (year)", 
             legend.width=1,
             slim=c(min(mij,na.rm=T), max(c(2,mij),na.rm=T)), 
             cex.lab=1.2,
             cex.main=1.2, 
             las=0.1,
             col=redc, 
             #main=paste("Average number of contacts per day",plot_title_extra), 
             #main=expression('m'['ij'] * .(plot_title_extra)), 
             main = bquote(paste('m'['ij']*' ', .(plot_title_extra))),
             xaxt="n", 
             yaxt="n")
  # set axis 
  plt_ticks <- seq(0,1,length=nrow(mij))
  axis(2, at=plt_ticks, labels = c(colnames(mij)),cex.axis=0.9,tick = FALSE,las=1)
  axis(1, at=plt_ticks, labels = c(colnames(mij)),cex.axis=0.9,tick = FALSE)
  
  # format results (rounding/scientific)
  if(any(mij>1,na.rm=T)){
    mij <- round(mij,digits=format_num_digits)
  } else{
    mij <- format(mij,digits = format_num_digits)
  }
  # get grid centers and add value
  e_grid <- expand.grid(plt_ticks,plt_ticks)
  text(e_grid, labels = mij)
}


