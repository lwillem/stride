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
#  Copyright 2024 Willem L, Kuylen E & Broeckhove J
############################################################################ #


if(0==1) # for debugging
{

  # load help functions
  source('./bin/rstride/Misc.R')
  
  #setwd('..')
  project_summary <- .rstride$load_project_summary(project_dir)
  exp_summary <- project_summary[2,]
  data_dir    <- './data'
  .rstride$plot_contacts(exp_summary,data_dir)
}

############################################################################# #
# INSPECT SOCIAL CONTACT PATTERNS                                          ####
############################################################################# # 

inspect_contact_data <- function(project_dir){
  
  # command line message
  smd_print('INSPECT SOCIAL CONTACT PATTERNS...')
  
  # load summary
  project_summary <- .rstride$load_project_summary(project_dir)
  
  # start slave nodes
  #smd_start_cluster()
  
  # analyse data
  i_exp <- 1
  # parallel processing issues to pass .rstride environment
  foreach(i_exp = 1:nrow(project_summary),
          .packages = 'simid.rtools',
          .export = '.rstride') %do% 
  {  
    # plot contacts
    .rstride$plot_contacts(project_dir,project_summary[i_exp,],'./data')
  }
  
  # end slave nodes
  smd_stop_cluster()
  
  # terminal message
  smd_print('INSPECTION OF SOCIAL CONTACTS PATTERNS COMPLETE')
}

############################################################################# #
# PLOT SOCIAL CONTACT MATRICES AND COUNTS                                  ####
############################################################################# #
# exp_summary <- project_summary[i_exp,]; data_dir <- './data'; bool_rds=FALSE;survey_start = exp_summary$start_date
.rstride$plot_contacts <- function(project_dir,exp_summary,data_dir)
{

  ##################### #
  ## GET DATA        ####
  ##################### #

  # load data
  data_cnt_all      <- .rstride$load_aggregated_output(project_dir,'data_contacts',exp_summary$exp_id)
  data_part_demo    <- .rstride$load_aggregated_output(project_dir,'data_participants',exp_summary$exp_id)
  data_health_all   <- .rstride$load_aggregated_output(project_dir,'data_health',exp_summary$exp_id)
  
  # select participant and contact data from from contact survey (and not infected seeds)
  data_part_demo  <- data_part_demo[data_part_demo$survey_type == "contacts",]
  data_cnt_all    <- data_cnt_all[data_cnt_all$local_id %in% data_part_demo$local_id,]
  data_health_all <- data_health_all[data_health_all$local_id %in% data_part_demo$local_id,]
  
  dim(data_part_demo)
  dim(data_health_all)
  
  # focus on demographic characteristics
  data_part_demo <- data_part_demo[,!grepl('is_',names(data_part_demo))]
  data_part_demo$sim_day <- NULL

  # join demographic and health data
  data_part_all <- merge(data_health_all,data_part_demo,by=c('local_id','exp_id'))

  ## reformat
  data_cnt_all$cnt_school    <- as.numeric(data_cnt_all$cnt_school)
  data_cnt_all$cnt_community_weekend <- as.numeric(data_cnt_all$cnt_community_weekend)
  data_cnt_all$cnt_community_weekday  <- as.numeric(data_cnt_all$cnt_community_weekday)
  data_cnt_all$sim_day       <- as.numeric(data_cnt_all$sim_day)
  data_cnt_all$cnt_prob      <- as.numeric(data_cnt_all$cnt_prob)
  data_cnt_all$part_sympt    <- as.numeric(data_cnt_all$part_sympt)
  data_cnt_all$cnt_sympt     <- as.numeric(data_cnt_all$cnt_sympt)
  
  summary(data_part_all)
  data_part_all$is_susceptible      <- as.numeric(data_part_all$is_susceptible)
  data_part_all$is_infected         <- as.numeric(data_part_all$is_infected)
  data_part_all$is_infectious       <- as.numeric(data_part_all$is_infectious)
  data_part_all$is_symptomatic      <- as.numeric(data_part_all$is_symptomatic)
  data_part_all$is_recovered        <- as.numeric(data_part_all$is_recovered)
  data_part_all$is_immune           <- as.numeric(data_part_all$is_immune)
  data_part_all$start_symptomatic   <- as.numeric(data_part_all$start_symptomatic)
  data_part_all$end_infectiousness  <- as.numeric(data_part_all$end_infectiousness)
  
  # if at least one data source is missing... stop
  if(nrow(data_cnt_all)==0 || nrow(data_part_all)==0) 
  {
    #smd_print("PARTICIPANT OR CONTACT DATA MISSING... STOP CONTACT ANALYSIS FOR",exp_summary$output_prefix)
    return(NULL) 
  }
    
  ## people without contacts
  dim(data_part_all)[1] - length(unique(data_cnt_all$local_id))
     
  ## employed and student population
  data_part_all$employed <- data_part_all$workplace_id != 0
  data_part_all$student  <- data_part_all$school_id != 0
  
  # get experiment tag
  exp_tag <- .rstride$create_exp_tag(exp_summary$exp_id)
  
  ## Socrates matrices ####
  
  # get results with default (minors and adults) or specific age groups    
  age_cat_breaks <- c(0,18,110)
  if('contact_survey_ages' %in% names(exp_summary)){
    age_cat_breaks <- as.numeric(unlist(strsplit(exp_summary$contact_survey_ages,',')))
  }
  plot_socrates_all(data_cnt_all   = data_cnt_all,
                    data_part_all  = data_part_all,
                    age_cat_breaks = age_cat_breaks,
                    project_dir    = project_dir,
                    exp_tag        = paste0(exp_tag,'_AG'),
                    survey_start   = exp_summary$start_date,
                    bool_rds       = TRUE)
  
  # # check rds
  # cnt_matrix_all <- readRDS(file.path(project_dir,'exp0001_AG_cnt_matrix.rds'))
  
  ## Other figures ####
  opt_day <- unique(data_cnt_all$sim_day)
  i_day <- 0
  for(i_day in opt_day){
    
    data_part_day <- data_part_all[data_part_all$sim_day == i_day,]
    data_cnt_day <- data_cnt_all[data_cnt_all$sim_day == i_day,]
    .rstride$plot_cnt_counts(data_part_day = data_part_day,
                             data_cnt_day  = data_cnt_day,
                             exp_tag       = exp_tag,
                             exp_summary   = exp_summary,
                             data_dir      = data_dir)
  }

  
  
  
} # end function

#################################  OTHER HELP FUNCTIONS  ################################ #


.rstride$plot_cnt_counts <- function(data_part_day,
                                     data_cnt_day,
                                     exp_tag,
                                     exp_summary,
                                     data_dir){
  
  ## SETTINGS 
  L <- max(c(80,data_part_day$part_age))
  num_days      <- length(unique(data_cnt_day$sim_day))
  sim_day       <- paste(unique(data_cnt_day$sim_day),collapse='-')
 
  # open pdf stream  
  .rstride$create_pdf(project_dir,paste0(exp_tag,'_cnt_patterns_d',sim_day),10,5)
  #par(mfrow=c(2,2))
  
  ## TOTAL
  mij_total  <- .rstride$plot_cnt_matrix(data_cnt_day,data_part_day,'total',L,num_days)
  
  ## HOUSEHOLD
  mij_hh     <- .rstride$plot_cnt_matrix(data_cnt_day[data_cnt_day$cnt_home==1,],data_part_day,'household',L,num_days)
  
  ## SCHOOL
  mij_school <- .rstride$plot_cnt_matrix(data_cnt_day[data_cnt_day$cnt_school==1,],data_part_day[data_part_day$student==T,],'school',L,num_days)
  
  ## WORKPLACE
  mij_workplace   <- .rstride$plot_cnt_matrix(data_cnt_day[data_cnt_day$cnt_workplace==1,],data_part_day[data_part_day$employed==T,],'workplace',L,num_days)
  
  ## WEEKEND COMMUNITY
  mij_community_weekend <- .rstride$plot_cnt_matrix(data_cnt_day[data_cnt_day$cnt_community_weekend==1,],data_part_day,'community_weekend',L,num_days)
  
  ## WEEKDAY COMMUNITY
  mij_community_weekday <- .rstride$plot_cnt_matrix(data_cnt_day[data_cnt_day$cnt_community_weekday==1,],data_part_day,'community_weekday',L,num_days)
  
  ## HOUSEHOLD CLUSTER
  .rstride$plot_cnt_matrix(data_cnt_day[data_cnt_day$cnt_household_cluster==1,],data_part_day,'household_cluster',L,num_days)
  
  #dev.off()
  
  ref_data_tag <- 'ref_fl2010'
  if(grepl('15touch',exp_summary$age_contact_matrix_file)){
    #ref_data_tag <- 'ref_fl2010_15touch'
    smd_print("NO REFERENCE 15_touch CONTACT DATA AVAIABLE",WARNING = TRUE)
  }
  
  # LOAD SURVEY DATA FROM FLANDERS AND FULLY CONNECTED HOUSEHOLDS
  survey_data <- xmlToList(file.path(data_dir,exp_summary$age_contact_matrix_file))
  names(survey_data)
  
  get_survey_data <- function(cluster_type,survey_data){
    survey_cluster     <- unlist(survey_data[[cluster_type]])
    flag_rate          <- grepl('contact.rate',names(survey_cluster))
    survey_mij_cluster <- matrix(as.numeric(survey_cluster[flag_rate]),nrow=sum(flag_rate))
    return(survey_mij_cluster)
  }
  
  survey_mij_hh         <- get_survey_data('household',survey_data)
  survey_mij_school     <- get_survey_data('school',survey_data)
  survey_mij_workplace  <- get_survey_data('workplace',survey_data)
  survey_mij_community  <- get_survey_data('community_weekday',survey_data)
  survey_mij_total      <- get_survey_data('regular_weekday',survey_data)
  
  survey_mij_school_weekend     <- survey_mij_school*0
  survey_mij_workplace_weekend  <- survey_mij_workplace*0
  survey_mij_community_weekend  <- get_survey_data('community_weekend',survey_data)
  survey_mij_total_weekend      <- get_survey_data('regular_weekend',survey_data)
  
  ## COMPARE
  par(mfrow=c(2,3))
  
  plot(rowSums(survey_mij_total),main='total',xlab='age',ylab='contacts',type='l',ylim=c(-0.1,40))
  lines(rowSums(survey_mij_total_weekend),main='total',xlab='age',ylab='contacts',type='l',lty=2)
  points(rowSums(mij_total,na.rm=T),col=2)
  legend('topright',c('week','weekend','model'),col=c(1,1,2),lty=c(1,2,0),pch=c(-1,-1,1),cex=0.8,title=ref_data_tag)
  
  plot(rowSums(survey_mij_hh),main='household',xlab='age',ylab='contacts',type='l',ylim=c(-0.1,5))
  points(rowSums(mij_hh,na.rm=T),col=2)
  legend('topright',c('week','weekend','model'),col=c(1,1,2),lty=c(1,2,0),pch=c(-1,-1,1),cex=0.8,title=ref_data_tag)
  
  plot(rowSums(survey_mij_school),main='school',xlab='age',ylab='contacts',type='l',ylim=c(-0.1,20))
  lines(rowSums(survey_mij_school_weekend),type='l',lty=2)
  points(rowSums(mij_school,na.rm=T),col=2)
  legend('topright',c('week','weekend','model'),col=c(1,1,2),lty=c(1,2,0),pch=c(-1,-1,1),cex=0.8,title=ref_data_tag)
  
  plot(rowSums(survey_mij_workplace),main='workplace',xlab='age',ylab='contacts',type='l',ylim=c(-0.1,20))
  lines(rowSums(survey_mij_workplace_weekend),type='l',lty=2)
  points(rowSums(mij_workplace,na.rm=T),col=2)
  legend('topright',c('week','weekend','model'),col=c(1,1,2),lty=c(1,2,0),pch=c(-1,-1,1),cex=0.8,title=ref_data_tag)
  
  plot(rowSums(survey_mij_community),main='weekend community',xlab='age',ylab='contacts',type='l',ylim=c(-0.1,25))
  lines(rowSums(survey_mij_community_weekend),type='l',lty=2)
  points(rowSums(mij_community_weekend,na.rm=T),col=2)
  legend('topright',c('week','weekend','model'),col=c(1,1,2),lty=c(1,2,0),pch=c(-1,-1,1),cex=0.8,title=ref_data_tag)
  
  plot(rowSums(survey_mij_community),main='weekday community',xlab='age',ylab='contacts',type='l',ylim=c(-0.1,25))
  lines(rowSums(survey_mij_community_weekend),type='l',lty=2)
  points(rowSums(mij_community_weekday,na.rm=T),col=2)
  legend('topright',c('week','weekend','model'),col=c(1,1,2),lty=c(1,2,0),pch=c(-1,-1,1),cex=0.8,title=ref_data_tag)
  par(mfrow=c(1,1))
  
  ## Transmission probability ####
  par(mfrow=c(2,2))
  cnt_location_opt <- c('cnt_home', 'cnt_school', 'cnt_workplace', 'cnt_community_weekend', 'cnt_community_weekday','cnt_household_cluster')
  for(i_cnt in cnt_location_opt){
    flag <- data_cnt_day[,i_cnt] == 1
    if(any(flag))
      boxplot(cnt_prob ~ part_age, data=data_cnt_day[flag,],
              main=paste(i_cnt, '[CNT]'),xlab='age',ylab='contact probability')
  }
  
  for(i_cnt in cnt_location_opt){
    flag <- data_cnt_day[,i_cnt] == 1
    if(any(flag))
      boxplot(trm_prob ~ part_age, data=data_cnt_day[flag,],
              main=paste(i_cnt, '[TRM]'),xlab='age',ylab='transmission probability')
  }
  
  dev.off() # close pdf stream
  
}


## RESHAPE DATA AND PLOT ####
#f_data_cnt = data_cnt_all;f_data_part=data_part_all;tag='total';L;num_days
.rstride$plot_cnt_matrix <- function(f_data_cnt,f_data_part,tag,L,num_days)
{
  
  # select participants
  data_cnt_flag <- f_data_cnt$local_id %in% f_data_part$local_id 
  
  # temporary max age
  L_temp <- max(f_data_cnt$part_age,L,f_data_cnt$cnt_age)+1
  
  # count contacts
  mij_tbl <- table(f_data_cnt$part_age,f_data_cnt$cnt_age)
  row_ind <- as.numeric(row.names(mij_tbl)) +1 # age 0 == index 1
  col_ind <- as.numeric(colnames(mij_tbl))  +1 # age 0 == index 1
  mij <- matrix(0,L_temp,L_temp)
  mij[row_ind,col_ind] <- mij_tbl
  
  # count participant per age
  if(sum(mij)==0)
  {
    row_ind <- 1:(L+1)
    col_ind <- 1:(L+1)
  }
  p_ages_tbl              <- table(f_data_part$part_age)
  
  row_ind <- as.numeric(names(p_ages_tbl)) +1
  p_ages <- matrix(0,L_temp,1)
  p_ages[row_ind] <- p_ages_tbl
  
  # remove age>L (age_column L+1)
  mij <- mij[1:(L+1),1:(L+1)]
  p_ages <- p_ages[1:(L+1)]
  
  # adjust for number of participants (if age is present)
  ages_present <- p_ages>0
  p_ages[ages_present]
  
  # contacts/participant
  for(j in 1:(L+1))
  {
    if(p_ages[j]>0)
    {
      mij[j,] <- mij[j,] / p_ages[j]
    }
  }
  
  # account for multiple days
  mij <- mij/num_days
  
  # set NA if participant age is not present in provided data
  #mij[,!ages_present] <- NA
  mij[!ages_present,] <- NA
  
  
  # plot matrix 
  g_matrix <- .rstride$plot_cnt_matrix_ggplot(mij,tag,FALSE)
  
  # plot number of contacts
  g_count  <- .rstride$plot_cnt_count_ggplot(f_data_cnt,f_data_part,L,num_days,tag)
  
  grid.arrange(g_matrix, g_count, ncol = 2)
  
  return(mij)
}

## PLOT CNT MATRIX ####
.rstride$plot_cnt_matrix_ggplot <- function(mij,title,bool_contour)
{
  ## remove small numbers
  mij[mij < quantile(mij,0.1,na.rm=T)] <- 0
  
  # Function to rescale z according to quantiles
  mij_ecdf <- ecdf(mij)
  
  # Covert matrix into data.frame for plotting with ggplot
  ggplot_data <- expand.grid(x = 0:(nrow(mij) - 1), y = 0:(nrow(mij) - 1))
  ggplot_data <- within(ggplot_data, {
    z <- as.vector(mij)
    z.rescaled <- mij_ecdf(mij)
  })
  
  
  z.breaks <- signif(unique(quantile(ggplot_data$z, prob = seq(from = 0, to = 1, length = 5),na.rm=T)), digits = 1)
  z.breaks.rescaled <- mij_ecdf(z.breaks)
  
  # Create the plot
  g <- ggplot(data = ggplot_data, mapping = aes(x = x, y = y, fill = z.rescaled, z = z.rescaled)) +
    geom_raster() +
    guides(colour=guide_legend("",order = 0),
           fill = guide_colourbar(order = 1))  +
    scale_fill_distiller(
      palette = "YlOrRd",
      breaks = z.breaks.rescaled,
      labels = z.breaks,
      name = 'Rate') +
    labs(x = "Age", y = "Age of contacts") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_fixed() +
    theme_bw() +
    ggtitle(title) +
    theme(legend.justification = c(1, 1),
          legend.position = 'right',
          legend.text     = element_text(size=18),
          legend.title    = element_text(size=18),
          axis.text       = element_text(size=20),
          axis.title      = element_text(size=20),
          plot.title      = element_text(size=40, face="bold",hjust = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()
    ) 
  
  # add missing ages?
  if(any(is.na(ggplot_data$z.rescaled))){
    g + geom_tile(data = subset(ggplot_data,  is.na(z.rescaled)), aes(colour = 'No data'),
                  linetype = 0, fill = "grey50")
  }
    # Add contour lines?
  if (bool_contour) {
    g + geom_contour(breaks = z.breaks, colour = "black", size = 0.2)
  } 
  
  # return
  return(g)
}

## PLOT CNT COUNT ####
.rstride$plot_cnt_count_ggplot <- function(f_data_cnt,f_data_part,L,num_days,title){
  
  if(nrow(f_data_cnt)==0){
    ggplot_data <-data.frame(local_id = -1,
                             part_age = f_data_part$part_age,
                             cnt_count = 0)
  } else{
    # Covert matrix into data.frame for plotting with ggplot
    ggplot_data        <- data.frame(table(f_data_cnt$local_id)/ num_days) 
    names(ggplot_data) <- c('local_id','cnt_count')
    ggplot_data <- merge(ggplot_data,f_data_part)
  }
  
  # remove oldest ages
  ggplot_data <- ggplot_data[ggplot_data$part_age<=L,]
  
  # average count per age
  cnt_age_mean <- aggregate(cnt_count~ part_age , data = ggplot_data ,mean) 
  
  # create plot
  g_plot <- ggplot(ggplot_data, aes(x=part_age, y=cnt_count)) + 
    aes(group = part_age) +
    geom_boxplot() +
    labs(x = "Age", y = "Count") +
    theme_bw() +
    ggtitle(title) +
    geom_line(data = cnt_age_mean, aes(x=part_age, y=cnt_count,group = 1),
              linewidth=2, colour="red") +
    theme(legend.justification = c(1, 1),
          legend.position = 'right',
          legend.text = element_text(size=18),
          legend.title = element_text(size=18),
          axis.text=element_text(size=20),
          axis.title=element_text(size=20),
          plot.title = element_text(size=40, face="bold",hjust = 0.5)
    )
  return(g_plot)
  
}


## SOCIAL CONTACT RATES OVER TIME ----
.rstride$contact_timeline <- function(project_dir)
{

  # load summary
  project_summary <- .rstride$load_project_summary(project_dir)
  
  # retrieve all variable model parameters
  input_opt_design     <- .rstride$get_variable_model_param(project_summary)
  input_opt_design
  
  # load data
  data_cnt      <- .rstride$load_aggregated_output(project_dir,'data_contacts')
  data_part     <- .rstride$load_aggregated_output(project_dir,'data_participants')

  ## reformat
  data_cnt$cnt_school              <- as.numeric(data_cnt$cnt_school)
  data_cnt$cnt_community_weekend  <- as.numeric(data_cnt$cnt_community_weekend)
  data_cnt$cnt_community_weekday  <- as.numeric(data_cnt$cnt_community_weekday)
  data_cnt$sim_day       <- as.numeric(data_cnt$sim_day)
  data_cnt$cnt_prob      <- as.numeric(data_cnt$cnt_prob)
  data_cnt$part_sympt    <- as.numeric(data_cnt$part_sympt)
  data_cnt$cnt_sympt     <- as.numeric(data_cnt$cnt_sympt)
  
  summary(data_part)
  data_part$is_susceptible      <- as.numeric(data_part$is_susceptible)
  data_part$is_infected         <- as.numeric(data_part$is_infected)
  data_part$is_infectious       <- as.numeric(data_part$is_infectious)
  data_part$is_recovered        <- as.numeric(data_part$is_recovered)
  data_part$is_immune           <- as.numeric(data_part$is_immune)
  data_part$start_symptomatic   <- as.numeric(data_part$start_symptomatic)
  data_part$end_infectiousness  <- as.numeric(data_part$end_infectiousness)
  
  
  
  data_cnt$sim_date <- data_cnt$sim_day + as.Date(unique(project_summary$start_date))
  
  # # specific summary
  # # summary_table                <- as.data.frame.matrix(table(data_transm[,colname_date],data_transm[,colname_value]))
  # summary_table                  <- dcast(data_cnt, formula('sim_date ~ exp_id'), value.var='ID', length)
  # names(summary_table)           <- c('sim_date',paste(prefix,names(summary_table)[-1],sep='_'))
  
  table(data_part$exp_id)
  
  
  ## COMMUNITY
  flag_cnt <- data_cnt$cnt_community_weekend == 1 | data_cnt$cnt_community_weekday == 1
  xx <- table(data_cnt$exp_id[flag_cnt],data_cnt$sim_date[flag_cnt])
  xx  / 8010
  
  cnt_dates <- as.Date(colnames(xx))
  com_cnt_pp <- xx / 8010
  
  flag_pre_lockdown <- cnt_dates < as.Date("2020-03-13")
  flag_lockdown     <- cnt_dates %in% as.Date("2020-03-14"):as.Date("2020-05-03")
  flag_exit_p1      <- cnt_dates %in% as.Date("2020-05-04"):as.Date("2020-05-09")
  flag_exit_p2      <- cnt_dates %in% as.Date("2020-05-10"):as.Date("2020-05-20")
  flag_exit_p3      <- cnt_dates %in% as.Date("2020-05-25"):as.Date("2020-05-27")
  
  
  
  mean(com_cnt_pp[,flag_pre_lockdown])
  mean(com_cnt_pp[,flag_lockdown])
  cnt_base <- mean(com_cnt_pp[,flag_lockdown])*(1/0.15)
  
  apply(com_cnt_pp[,flag_exit_p1],1,mean) / cnt_base
  apply(com_cnt_pp[,flag_exit_p2],1,mean) / cnt_base
  apply(com_cnt_pp[,flag_exit_p3],1,mean) / cnt_base
  
  
  plot(as.Date(colnames(xx)),
       xx[1,]/ 8010 ,ylim=range(xx/ 8010))
  points(as.Date(colnames(xx)),
       xx[2,]/ 8010 )
  points(as.Date(colnames(xx)),
         xx[3,]/ 8010)
  points(as.Date(colnames(xx)),
         xx[4,]/ 8010)
  
  ## COMBINED
  flag_cnt <- data_cnt$cnt_community_weekend == 1 | data_cnt$cnt_community_weekday == 1 | data_cnt$cnt_household_cluster == 1
  yy <- table(data_cnt$exp_id[flag_cnt],data_cnt$sim_date[flag_cnt])
  plot(as.Date(colnames(yy)),
         yy[3,]/ 8010,
         col=3) 
  points(as.Date(colnames(yy)),
         yy[2,]/ 8010,
         col=2) 
  
  points(as.Date(colnames(yy)),
         yy[1,]/ 8010,
         col=1) 
  
 ## HH BUBBLE
  flag_cnt <- data_cnt$cnt_household_cluster == 1
  zz <- table(data_cnt$exp_id[flag_cnt],data_cnt$sim_date[flag_cnt])
    
}


