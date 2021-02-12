############################################################################ #
#  This file is part of the Stride software. 
#
#  Copyright 2021, Willem L
############################################################################ #
#
# TO CREATE CALENDAR FILE(S) FOR 2019-2021
#
# CONTAINING:
# 1. Public and school holidays
# 2. Contact reductions (covid19)
#       * pre-, primairy and secondary school
#       * workplace
#       * community
#       * household clusters
# 3. Imported cases (covid19)
# 4. Contact tracing (covid19)
# 5. Universal testing (covid19)
#
############################################################################ #

# debug
if(0==1){
  source('bin/rstride/rStride.R')
  cnt_other_exit_delay <- 21
  show_plots = TRUE
  
  create_calendar_file(show_plots=TRUE)
}


# create calendar files if they do not exist, else re-use them
create_calendar_file <- function(file_name_tag='2020_2021',show_plots = FALSE)
{
  
  filename_calendar_full <- paste('sim_output/calendar_belgium',file_name_tag,'covid19.csv',sep='_')

  ########################################### #
  ## INITIATE DATA                       ####
  ########################################### #
  
  date_start <- as.Date("2019-01-01")
  date_end   <- as.Date("2021-12-31")
  
  # default value in C++ CALENDAR vectors ==>> 0
  
  ########################################### #
  ## 1.a Public holidays                 ####
  ########################################### #
  
  data.table(category = "general",
             date     = as.Date(c(
               '2019-01-01','2019-04-22','2019-05-01','2019-05-30','2019-06-10', # 2019
               '2019-07-21','2019-08-15','2019-11-01','2019-11-11','2019-12-25',
               
               '2020-01-01','2020-04-13','2020-05-01','2020-05-21','2020-06-01', # 2020
               '2020-07-21','2020-08-15','2020-11-01','2020-11-11','2020-12-25',
               
               '2021-01-01','2021-04-05','2021-05-01','2021-05-13','2021-06-24', # 2021
               '2021-07-21','2021-08-15','2021-11-01','2021-11-11','2021-12-25')),
             value    = 1,
             type = 'boolean',
             age = NA_integer_,
             stringsAsFactors = F
            ) -> d_calendar_holiday
  
  summary(d_calendar_holiday)
  
  ########################################### #
  ## 1.b School holidays                 ####
  ########################################### #
  
  # start from blanc calendar
  data.table(category = "schools_closed",
             date     = seq(as.Date(date_start),as.Date(date_end),1),
             value    = 0.0,
             type = 'double',
             age = NA_integer_,
             stringsAsFactors = F
  ) -> d_school_holidays
  d_college_holidays <- copy(d_school_holidays)
  
  # include K12 school holidays
  d_school_holidays[date %in% c(seq(as.Date('2019-01-01'),as.Date('2019-01-06'),1), # 2019
                                seq(as.Date('2019-03-04'),as.Date('2019-03-10'),1),
                                seq(as.Date('2019-04-08'),as.Date('2019-04-22'),1),
                                seq(as.Date('2019-07-01'),as.Date('2019-08-31'),1),
                                seq(as.Date('2019-10-28'),as.Date('2019-11-03'),1),
                                seq(as.Date('2019-12-23'),as.Date('2019-12-31'),1),
                                
                                seq(as.Date('2020-01-01'),as.Date('2020-01-05'),1), # 2020
                                seq(as.Date('2020-02-24'),as.Date('2020-02-29'),1),
                                seq(as.Date('2020-04-06'),as.Date('2020-04-19'),1),
                                seq(as.Date('2020-07-01'),as.Date('2020-08-31'),1),
                                seq(as.Date('2020-11-02'),as.Date('2020-11-08'),1),
                                seq(as.Date('2020-12-21'),as.Date('2020-12-31'),1),
                                
                                seq(as.Date('2021-01-01'),as.Date('2021-01-03'),1), # 2021
                                seq(as.Date('2021-02-15'),as.Date('2021-02-21'),1),
                                seq(as.Date('2021-04-05'),as.Date('2021-04-18'),1),
                                seq(as.Date('2021-07-01'),as.Date('2021-08-31'),1),
                                seq(as.Date('2021-11-01'),as.Date('2021-11-07'),1),
                                seq(as.Date('2021-12-27'),as.Date('2021-12-31'),1)),
             value    := 1.0]
  
  # add college holidays
  d_college_holidays[date %in% c(seq(as.Date('2019-01-01'),as.Date('2019-01-06'),1), # 2019
                          seq(as.Date('2019-03-04'),as.Date('2019-03-10'),1),
                          seq(as.Date('2019-04-08'),as.Date('2019-04-22'),1),
                          seq(as.Date('2019-07-01'),as.Date('2019-09-22'),1), # summer break untill September, 22
                          #seq(as.Date('2019-10-28'),as.Date('2019-11-03'),1), # no fall break
                          seq(as.Date('2019-12-23'),as.Date('2019-12-31'),1),
                          
                          seq(as.Date('2020-01-01'),as.Date('2020-01-05'),1), # 2020
                          seq(as.Date('2020-02-24'),as.Date('2020-02-29'),1),
                          seq(as.Date('2020-04-06'),as.Date('2020-04-19'),1),
                          seq(as.Date('2020-07-01'),as.Date('2020-09-20'),1),# summer break untill September, 20
                          #seq(as.Date('2020-11-02'),as.Date('2020-11-08'),1), # no fall break
                          seq(as.Date('2020-12-21'),as.Date('2020-12-31'),1),
                          
                          seq(as.Date('2021-01-01'),as.Date('2021-01-03'),1), # 2021
                          seq(as.Date('2021-02-15'),as.Date('2021-02-21'),1),
                          seq(as.Date('2021-04-05'),as.Date('2021-04-18'),1),
                          seq(as.Date('2021-07-01'),as.Date('2021-09-19'),1), # summer break untill September, 19
                          #seq(as.Date('2021-11-01'),as.Date('2021-11-07'),1), # no fall break
                          seq(as.Date('2021-12-27'),as.Date('2021-12-31'),1)),
             value    := 1.0]
  
  #K12 school
  tmp_school_holidays <- copy(d_school_holidays)
  tmp_school_holidays[,category:='schools_closed']
  for(i_age in 0:17){
    d_calendar_holiday <- rbind(d_calendar_holiday,copy(tmp_school_holidays[,age:=i_age]))
  }
  
  # College
  tmp_college_holidays <- copy(d_college_holidays)
  tmp_college_holidays[,category:='schools_closed']
  for(i_age in 18:25){
    d_calendar_holiday <- rbind(d_calendar_holiday,copy(tmp_college_holidays[,age:=i_age]))
  }
  
  
  ########################################################### #
  ##  2a. Contact reductions: school closures              ####
  ########################################################### #
  #       * (pre-, primary and secondary school)
  
  # set default school closure
  # school_dates_non_holiday <- date_all[!date_all %in% dcal_school_closure]
  # data.table(category = "schools_closed",
  #            # date     = seq(as.Date('2020-03-14'),as.Date('2020-06-30'),1),
  #            date     = school_dates_non_holiday,
  #            value    = 0.0,
  #            type = 'double',
  #            age = NA_integer_,
  #            stringsAsFactors = F
  # ) -> d_school_closure
  # 
  # d_school_closure[date %in% seq(as.Date('2020-03-14'),as.Date('2020-06-30'),1),value := 1.0]
  # 
  # tmp_school_closure <- copy(d_school_closure)
  # tmp_school_closure[,category:='schools_closed']
  
  d_calendar_holiday[date %in% seq(as.Date('2020-03-14'),as.Date('2020-06-30'),1), value := 1.0] 

  # set eligible dates for school reopening in May-June 2020
  d_school_reopening <- seq(as.Date('2020-05-18'),as.Date('2020-06-30'),1)
  d_school_reopening_wday <- as.POSIXlt(d_school_reopening)$wday
  
  # preschool (reopens 4d/week)
  d_school_reopening_4d <- d_school_reopening[d_school_reopening_wday %in% 1:4]
  d_calendar_holiday[date %in% d_school_reopening_4d & age %in% c(0,1,2,6,7),value:=0.5]
  
  # primary school (reopens 2d/week)
  d_school_reopening_2d <- d_school_reopening[d_school_reopening_wday %in% 4:5]
  d_school_reopening_2d[1:2] <- d_school_reopening_2d[1:2] - 2 # fix for holidays Thu-Fri in May
  d_calendar_holiday[date %in% d_school_reopening_2d & age %in% c(11),value:=0.5]
  
  
  #secondary school (reopens 1d week)
  d_school_reopening_1d <- d_school_reopening[d_school_reopening_wday %in% 3]
  d_calendar_holiday[date %in% d_school_reopening_1d & age %in% c(17),value:=0.5]
  
  ########################################### #
  ##  2b. Contact reductions: other        ####
  ########################################### #
  #       * workplace
  #       * community
  #       * household clusters
  
  # workplace distancing
  data.table(category = "workplace_distancing",
             #date     = seq(as.Date('2020-03-14'),as.Date('2020-05-03'),1),
             date     = seq(as.Date(date_start),as.Date(date_end),1),
             value    = 0.0,
             type     = 'double',
             age = NA_integer_,
             stringsAsFactors = F
  ) -> dcal_workplace_distancing
  dcal_workplace_distancing[date %in% seq(as.Date('2020-03-14'),as.Date('2020-05-03'),1),value := 1.0]
  
  
  # community distancing
  data.table(category = "community_distancing",
             #date     = seq(as.Date('2020-03-14'),as.Date('2020-05-24'),1),
             date     = seq(as.Date(date_start),as.Date(date_end),1),
             value    = 0.0,
             type = 'double',
             age = NA_integer_,
             stringsAsFactors = F
  ) -> dcal_community_distancing
  dcal_community_distancing[date %in% seq(as.Date('2020-03-14'),as.Date('2020-05-24'),1),value := 1.0]
  
  # household clustering
  data.table(category = "household_clustering",
             date     = seq(as.Date('2020-05-11'),as.Date('2020-08-31'),1),
             value    = 1,
             type = 'boolean',
             age = NA_integer_,
             stringsAsFactors = F
  ) -> dcal_household_clustering
  
  ########################################### #
  ## 3. Imported cases                 ####
  ########################################### #
  
  data.table(category = "imported_cases",
             date     = seq(as.Date('2020-07-01'),as.Date('2020-08-31'),1),
             value    = 1,
             type = 'boolean',
             age = NA_integer_,
             stringsAsFactors = F
  ) -> dcal_imported_cases
  
  ########################################### #
  ##  4. Contact tracing                 ####
  ########################################### #
  
  data.table(category = "contact_tracing",
             date     = seq(as.Date('2020-05-11'),as.Date('2020-08-31'),1),
             value    = 1,
             type = 'boolean',
             age = NA_integer_,
             stringsAsFactors = F
  ) -> dcal_contact_tracing
  
  ########################################### #
  ## 5. Universal testing                 ####
  ########################################### #

  data.table(category = "universal_testing",
             date     = seq(as.Date('2020-05-11'),as.Date('2020-08-31'),1),
             value    = 1,
             type = 'boolean',
             age = NA_integer_,
             stringsAsFactors = F
  ) -> dcal_universal_testing
  
  
  ########################################### #
  ## MERGE HOLIDAYS & OTHER CALENDAR ITEMS ####
  ########################################### #
  
  # get 'dcal_*' variables
  opt_other <- ls(pattern='dcal_')
  
  # combine all 'dcal_*' variable
  d_calendar_all <- foreach(i_other  = opt_other,
                            .init    = d_calendar_holiday,
                            .combine = 'rbind') %do% {
                              get(i_other)
                            } 
  
  # select range
  d_calendar_all <- d_calendar_all[date >= date_start & date <= date_end,]
  range(d_calendar_all$date)
  
  ########################################### #
  ## EXPLORE DATA                         ####
  ########################################### #
  
  # open pdf stream
  pdf(file=gsub('.csv','.pdf',filename_calendar_full),6,6)
  
  plot_calendar(d_calendar_all,show_plots)
  
  # close pdf stream
  dev.off()
  
  ########################################### #
  ## SAVE AS CSV	 	         ####
  ########################################### #
  
  # # format date
  # d_calendar_all[,date:=format(date,'%Y-%m-%d')]
  # format(d_calendar_all$date,'%Y-%m-%d')
  
  # save as csv (all calendar info)
  write.table(d_calendar_all,
              file = filename_calendar_full,sep=',',row.names=F,quote=F)
  
  unique(d_calendar_all$category)
  
  return(filename_calendar_full)
}

# note: variable "b_school_repopening" is not used anymore... but still here for backward compatibility
plot_calendar <- function(dt_calendar, show_plots = TRUE, b_school_repopening=TRUE){
  if(show_plots){
    category_opt <- unique(dt_calendar$category)
    par(mfrow=c(4,2))
    
    # check if dt_calendar is data.table
    if(!is.data.table(dt_calendar)){
      dt_calendar <- data.table(dt_calendar)
    }
    
    # make sure that "date" is in date format
    dt_calendar$date <- as.Date(dt_calendar$date)
    
    x_lim <- range(dt_calendar$date)
    x_lab_year <- paste(unique(year(dt_calendar$date)),sep='-')
    i_cat <- category_opt[2]
    for(i_cat in category_opt){
      plot(x   = dt_calendar[category == i_cat,date],
           y   = dt_calendar[category == i_cat,value],
           xlim = x_lim,
           ylim = range(0,1,dt_calendar$value),
           col  = 1,
           #type='l',
           pch  = 15,
           #lwd=2,
           main = i_cat,
           bty='n',
           xlab = x_lab_year,
           ylab = unique(dt_calendar[,type]),
           xaxt = 'n'
      )
      add_x_axis(x_lim)
    }
    
    if("schools_closed" %in% dt_calendar$category){
      i_cat <- "schools_closed"
      plot(x   = dt_calendar[category == i_cat & value == 1,date],
           y   = dt_calendar[category == i_cat & value == 1,age],
           xlim = x_lim,
           ylim = range(0,1,dt_calendar$age,na.rm=T),
           col  = 1,
           pch  = 15,
           main = i_cat,
           bty='n',
           xlab = x_lab_year,
           ylab = 'age',
           xaxt = 'n'
      )
      points(x   = dt_calendar[category == i_cat & !value %in% 0:1,date],
           y   = dt_calendar[category == i_cat & !value %in% 0:1,age],
           col  = 13,
           pch  = 15
      )
      add_x_axis(x_lim)
    }
  }
}

# db_cat <- "workplace_distancing";date_start<-"2020-03-14";date_end <- "2020-05-02";db_value<-0.75;file_name<-dcal_file
# db_cat <- "workplace_distancing";date_selection <- wp_dist$x; db_value <- wp_dist$y; file_name <- dcal_file
# db_update <- data.frame(c('2020-02-10',0),
#                         c('2020-03-13',0),
#                         c('2020-03-19',0.85),
#                         c('2020-05-02',0.85),
# c('2020-05-03',0.75)); db_cat <- "workplace_distancing";db_age = 'NA'; file_name <- "sim_output/calendar_belgium_wp_fitting_covid19.csv"
# db_category =  "workplace_distancing";db_update = dcal_wp_distancing; file_name = dcal_file
adjust_calendar_file <- function(db_category, db_update, file_name, db_age = 'NA', show_plots=FALSE){
  
  # read calendar file
  d_calendar_all <- data.table(read.table(file=file_name,sep=',',header=T))
  
  # make sure the value is stored as double
  d_calendar_all$value <- as.double(d_calendar_all$value)
  
  # adjust age == NA
  d_calendar_all[, age_char := as.character(age)]
  d_calendar_all[is.na(age_char), age_char := 'NA']
  
  # check category
  if(!db_category %in% unique(d_calendar_all$category)){
    smd_print("CALENDAR CATEGORY UNKNOWN => STOP CALENDAR ADJUSTMENT")
    smd_print("CALENDAR CATEGORY OPTIONS:", paste0(unique(d_calendar_all$category),collapse=', '))
    
    return(NA)
  }
  
  # extrapolate given dates and values
  db_update  <- data.frame(t(db_update))
  date_out   <- seq(min(as.Date(db_update[,1])),max(as.Date(db_update[,1])),1) 
  db_update  <- approx(x=as.Date(db_update[,1]),
                      y=db_update[,2],
                      xout = date_out)
  names(db_update) <- c('date','value')

  # update values for all given date and ages
  for(i_db_age in as.character(db_age)){
    d_calendar_all[as.character(date) %in% as.character(db_update$date) &
                      category == db_category &
                      age_char == i_db_age, 
                   value := db_update$value ]
  }

  # explore
  plot_calendar(d_calendar_all,show_plots = show_plots)
  
  # save as csv 
  write.table(d_calendar_all,
              file = file_name,sep=',',row.names=F,quote=F)
  
}


if(0==1){ # debug----
  
  dcal_file <- create_calendar_file(file_name_tag = 'wp_fitting',show_plots = T)
  dcal_wp_distancing <- data.frame(c('2020-03-13',0),
                                   c('2020-03-19',0.85),
                                   c('2020-05-02',0.85),
                                   c('2020-05-03',0.5))
  
  adjust_calendar_file(db_category =  "workplace_distancing",db_update = dcal_wp_distancing, file_name = dcal_file)
  
  
}

