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
plot_socrates_all <- function(data_cnt_all,
                              data_part_all,
                              age_cat_breaks,
                              project_dir,
                              exp_tag,
                              survey_start,
                              bool_rds=FALSE){
  
  # open pdf stream
  .rstride$create_pdf(project_dir,paste0(exp_tag,'_cnt_matrix_all'),12,6)
  
  opt_day <- unique(data_cnt_all$sim_day)
  
  # option to initiate list to store contact data
  if(bool_rds){
    mij_summary <- list()
  }
  
  i_day <- 43
  for(i_day in opt_day){
    data_cnt_day <- data_cnt_all[data_cnt_all$sim_day == i_day,]
    mij_all      <- plot_socrates_location(data_cnt_day,data_part_all,age_cat_breaks,as.Date(survey_start) + i_day)
 
    # select symptomatic participants and their contacts on 'i_day' 
    # note: symptomatic people are only identified if they have contacts on day i_day (potential bias!)
    data_cnt_day    <- data_cnt_all[data_cnt_all$sim_day == i_day & data_cnt_all$part_sympt == 1,]
    data_part_sympt <- data_part_all[data_part_all$local_id %in% data_cnt_day$local_id,]
    mij_sympt       <- plot_socrates_location(data_cnt_day,data_part_sympt,age_cat_breaks,as.Date(survey_start) + i_day,title_add='SYMPT')
  
    # select non-symptomatic participants and their contacts on 'i_day' 
    data_cnt_day        <- data_cnt_all[data_cnt_all$sim_day == i_day & data_cnt_all$part_sympt == 0,]
    data_part_non_sympt <- data_part_all[!data_part_all$local_id %in% data_part_sympt$local_id,]
    mij_non_sympt       <- plot_socrates_location(data_cnt_day,data_part_non_sympt,age_cat_breaks,as.Date(survey_start) + i_day,title_add='NON-SYMPT')

    if(bool_rds){
      mij_summary[[paste0('day',i_day)]] <- list(mij_all       = mij_all,
                                                 mij_sympt     = mij_sympt,
                                                 mij_non_sympt = mij_non_sympt,
                                                 date          = as.Date(survey_start) + i_day,
                                                 exp_tag       = exp_tag)
    } 
    
  } # end for-loop opt_days
  dev.off()
  
  # option to store the contact matrices as rds file
  if(bool_rds){
    saveRDS(mij_summary,file=smd_file_path(project_dir,paste0(exp_tag,'_cnt_matrix.rds')))
  }
}

 # data_cnt <- data_cnt_day; data_part <- data_part_all;survey_day <- as.Date(survey_start) + i_day;title_add=''
plot_socrates_location <- function(data_cnt,data_part,age_cat_breaks,survey_day,title_add=''){
  
  par(mfrow=c(2,3))
  
  ## TOTAL
  mij_total <- plot_contact_matrix_socrates(data_cnt,data_part,paste(title_add,'total'),age_cat_breaks)
  
  ## HOUSEHOLD
  mij_household <- plot_contact_matrix_socrates(data_cnt[data_cnt$cnt_home==1,],data_part,paste(title_add,'@household'),age_cat_breaks)
  
  ## SCHOOL
  mij_school             <- plot_contact_matrix_socrates(data_cnt[data_cnt$cnt_school==1,],data_part,paste(title_add,'@school'),age_cat_breaks)
  mij_school_conditional <- plot_contact_matrix_socrates(data_cnt[data_cnt$cnt_school==1,],data_part[data_part$student==T,],paste(title_add,'@school (conditional)'),age_cat_breaks,bool_plot = FALSE)
  
  ## WORK
  mij_workplace             <- plot_contact_matrix_socrates(data_cnt[data_cnt$cnt_work==1,],data_part,paste(title_add,'@work'),age_cat_breaks)
  mij_workplace_conditional <- plot_contact_matrix_socrates(data_cnt[data_cnt$cnt_work==1,],data_part[data_part$employed==T,],paste(title_add,'@work (conditional)'),age_cat_breaks,bool_plot = FALSE)
  
  ## PRIMARY COMMUNITY
  mij_community_weekend <- plot_contact_matrix_socrates(data_cnt[data_cnt$cnt_community_weekend==1,],data_part,paste(title_add,'@weekend community'),age_cat_breaks)
  
  ## SECONDARY COMMUNITY
  mij_community_weekday <- plot_contact_matrix_socrates(data_cnt[data_cnt$cnt_community_weekday==1,],data_part,paste(title_add,'@week community'),age_cat_breaks)
  
  ## HOUSEHOLD CLUSTER
  mij_household_cluster <- plot_contact_matrix_socrates(data_cnt[data_cnt$cnt_household_cluster==1,],data_part,paste(title_add,'@household cluster'),age_cat_breaks)
  
  plot(0,0,col=0,axes=F,xlab='',ylab='')
  # survey_day <- as.Date(exp_summary$start_date) + unique(data_cnt$sim_day)
  survey_day <- as.Date(survey_day)
  Sys.setlocale("LC_TIME", 'en_GB.UTF-8')
  text(0,0,paste(title_add,format(survey_day,'\n%A'),survey_day),pos=3)
  text(0,0,paste('Contacts when symptomatic:',sum(data_cnt$part_sympt),
                  '\nNumber of participants:',nrow(data_part)),pos=1)
  
  return(list(mij_total                 = mij_total$matrix,
              mij_household             = mij_household$matrix,
              mij_school                = mij_school$matrix,
              mij_school_conditional    = mij_school_conditional$matrix,
              mij_workplace             = mij_workplace$matrix,
              mij_workplace_conditional = mij_workplace_conditional$matrix,
              mij_community_weekend     = mij_community_weekend$matrix,
              mij_community_weekday     = mij_community_weekday$matrix,
              mij_household_cluster     = mij_household_cluster$matrix,
              participants              = mij_total$participants))
}

# data_cnt <- data_cnt[data_cnt$cnt_work==1,]; data_part <- data_part[data_part$student==T,]
plot_contact_matrix_socrates <- function(data_cnt,data_part,figure_title,age_cat_breaks,bool_plot = TRUE){
  
   # select participant id and age
   data_part_sel <- data_part[,c('local_id','part_age')]
 
   # make sure the participant data is not empty and contains the maximum age
   # this dummy participant has no contacts, so has no impact on the final results
   if(nrow(data_part_sel)==0 || max(data_part_sel$part_age < max(age_cat_breaks))){
     data_part_sel <- rbind(data_part_sel,c(local_id=NA,part_age=max(age_cat_breaks)))
     colnames(data_part_sel) <- c('local_id','part_age')
   }
      
  # get socialmixr 'participants' object
  db_participants   <- data.frame(part_id     = data_part_sel$local_id,
                                  part_age    = data_part_sel$part_age,
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
                                              "cnt_community_weekend","cnt_community_weekday","part_sympt","cnt_sympt" )])
  
  # get socialmixr 'survey' object
  survey_rstride <- survey(participants = db_participants,
                           contacts     = db_contacts)

  # get matrix
  suppressWarnings(
  cnt_matrix <- contact_matrix(survey_rstride,age.limits = age_cat_breaks)  
  )
  
  # account for NA
  cnt_matrix$matrix[is.na(cnt_matrix$matrix)] <- 0
  
  # plot matrix
  if(bool_plot && any(cnt_matrix$matrix>0)){
    plot_cnt_matrix(cnt_matrix$matrix,figure_title)
  }
  
  # return matrix
  return(cnt_matrix)
}

#mij <- contact_matrix(polymod, countries = "United Kingdom", age.limits = c(0, 1, 5, 15))$matrix
#mij <- cnt_matrix$matrix
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
  if(any(mij>1e-2,na.rm=T)){
    mij_labels <- round(mij,digits=format_num_digits)
    cex_labels  <- 1
  } else{
    mij_labels <- format(mij,digits = format_num_digits)
    cex_labels <- 0.5
  }
  # get grid centres and add value
  e_grid <- expand.grid(plt_ticks,plt_ticks)
  text(e_grid, labels = mij_labels,cex = cex_labels)
}


