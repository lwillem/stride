############################################################################ #
#  This file is part of the Stride software. 
#
#  Copyright 2024
############################################################################ #
#
# TO CREATE CALENDAR FILE(S) FOR 2019-2021
#
# CONTAINING:
# 1. Public and school holidays
# 2. Contact reductions (covid19)
#       * pre-, primary and secondary school
#       * workplace
#       * community
#       * household clusters
# 3. Imported cases (covid19)
# 4. Contact tracing (covid19)
#
############################################################################ #

# debug
if(0==1){
  source('bin/rstride/rStride.R')
  show_plots = TRUE
  create_calendar_file(show_plots=TRUE)
}


# create calendar files if they do not exist, else re-use them
create_calendar_file <- function(file_name_tag = '2023_2026',
                                 date_end = "2026-12-31",
                                 show_plots = FALSE,
                                 file_name = NA)
{
  
  filename_calendar_full <- ifelse(is.na(file_name),
                                   paste('sim_output/calendar_dane',file_name_tag,'measles.csv',sep='_'),
                                   file_name)

  ########################################### #
  ## INITIATE DATA                       ####
  ########################################### #
  
  date_start <- as.Date("2023-01-01")

  # default value in C++ CALENDAR vectors ==>> 0
  
  ########################################### #
  ## 1.a Public holidays                 ####
  ########################################### #
  
  data.table(category = "general",
             date     = as.Date(c(
               # 2023
               '2023-01-02','2023-01-16','2023-02-20','2023-05-29','2023-06-19', # Jan 1 fell on Sunday
               '2023-07-04','2023-09-04','2023-10-09','2023-11-10','2023-11-23','2023-12-25', # Nov 11 fell on Saturday
               # 2024
               '2024-01-01','2024-01-15','2024-02-19','2024-05-27','2024-06-19', 
               '2024-07-04','2024-09-02','2024-10-14','2024-11-11','2024-11-28','2024-12-25',
               # 2025
               '2025-01-01','2025-01-20','2025-02-17','2025-05-26','2025-06-19', 
               '2025-07-04','2025-09-01','2025-10-13','2025-11-11','2025-11-27','2025-12-25',
               # 2026
               '2026-01-01','2026-01-19','2026-02-16','2026-05-25','2026-06-19',
               '2026-07-03','2026-09-07','2026-10-12','2026-11-11','2026-11-26','2026-12-25'
             )),
             value    = 1,
             type = 'boolean',
             age = NA_integer_,
             stringsAsFactors = F
            ) -> d_calendar_holiday
  
  summary(d_calendar_holiday)
  
  ########################################### #
  ## 1.b School holidays                 ####
  ########################################### #
  
  # start from blank calendar
  data.table(category = "schools_closed",
             date     = seq(as.Date(date_start),as.Date(date_end),1),
             value    = 0.0,
             type = 'double',
             age = NA_integer_,
             stringsAsFactors = F
  ) -> d_school_holidays
  d_college_holidays <- copy(d_school_holidays)
  
  # include school holidays
  d_school_holidays[date %in% c(
    # --- 2022-23 school year (source: official MMSD PDF) ---
    # seq(as.Date('2022-09-05'), as.Date('2022-09-05'), 1),  # Labor Day
    # seq(as.Date('2022-11-23'), as.Date('2022-11-25'), 1),  # Fall/Thanksgiving Break
    # seq(as.Date('2022-12-21'), as.Date('2023-01-03'), 1),  # Winter Break
    seq(as.Date('2023-01-01'), as.Date('2023-01-03'), 1),  # Winter Break
    # seq(as.Date('2023-01-16'), as.Date('2023-01-16'), 1),  # MLK Day
    seq(as.Date('2023-03-27'), as.Date('2023-03-31'), 1),  # Spring Break
    # seq(as.Date('2023-05-29'), as.Date('2023-05-29'), 1),  # Memorial Day
    seq(as.Date('2023-06-09'), as.Date('2023-08-31'), 1),  # Summer Break
    
    # --- 2023-24 school year (EXTRAPOLATED - no official MMSD calendar available;
    #     dates below based on MMSD's typical pattern + third-party aggregator data) ---
    # seq(as.Date('2023-09-04'), as.Date('2023-09-04'), 1),  # Labor Day
    seq(as.Date('2023-11-22'), as.Date('2023-11-24'), 1),  # Fall/Thanksgiving Break (extrapolated)
    seq(as.Date('2023-12-20'), as.Date('2024-01-02'), 1),  # Winter Break (extrapolated)
    # seq(as.Date('2024-01-15'), as.Date('2024-01-15'), 1),  # MLK Day
    seq(as.Date('2024-03-25'), as.Date('2024-03-29'), 1),  # Spring Break (extrapolated)
    # seq(as.Date('2024-05-27'), as.Date('2024-05-27'), 1),  # Memorial Day
    seq(as.Date('2024-06-07'), as.Date('2024-08-31'), 1),  # Summer Break (extrapolated)
    
    # --- 2024-25 school year (source: official MMSD calendar page) ---
    # seq(as.Date('2024-09-02'), as.Date('2024-09-02'), 1),  # Labor Day
    # seq(as.Date('2024-10-03'), as.Date('2024-10-04'), 1),  # Mid-Fall Break - MADISON SPECIFIC
    seq(as.Date('2024-11-27'), as.Date('2024-11-29'), 1),  # Fall/Thanksgiving Break
    seq(as.Date('2024-12-23'), as.Date('2025-01-03'), 1),  # Winter Break
    # seq(as.Date('2025-01-20'), as.Date('2025-01-20'), 1),  # MLK Day
    seq(as.Date('2025-03-24'), as.Date('2025-03-28'), 1),  # Spring Break
    # seq(as.Date('2025-05-26'), as.Date('2025-05-26'), 1),  # Memorial Day
    seq(as.Date('2025-06-12'), as.Date('2025-08-31'), 1),  # Summer Break
    
    # --- 2025-26 school year (source: official MMSD calendar page) ---
    # seq(as.Date('2025-09-01'), as.Date('2025-09-01'), 1),  # Labor Day
    # seq(as.Date('2025-10-16'), as.Date('2025-10-17'), 1),  # Mid-Fall Break - MADISON SPECIFIC
    seq(as.Date('2025-11-26'), as.Date('2025-11-28'), 1),  # Fall/Thanksgiving Break
    seq(as.Date('2025-12-22'), as.Date('2026-01-02'), 1),  # Winter Break
    # seq(as.Date('2026-01-19'), as.Date('2026-01-19'), 1),  # MLK Day
    seq(as.Date('2026-03-23'), as.Date('2026-03-27'), 1),  # Spring Break
    # seq(as.Date('2026-05-25'), as.Date('2026-05-25'), 1),  # Memorial Day
    seq(as.Date('2026-06-11'), as.Date('2026-08-31'), 1),  # Summer Break
    
    # --- 2026-27 school year (source: official MMSD calendar page) ---
    # seq(as.Date('2026-09-07'), as.Date('2026-09-07'), 1),  # Labor Day
    # seq(as.Date('2026-10-15'), as.Date('2026-10-16'), 1),  # Mid-Fall Break - MADISON SPECIFIC
    seq(as.Date('2026-11-25'), as.Date('2026-11-27'), 1),  # Fall/Thanksgiving Break
    seq(as.Date('2026-12-21'), as.Date('2026-12-31'), 1)  # Winter Break
    # seq(as.Date('2027-01-18'), as.Date('2027-01-18'), 1),  # MLK Day
    # seq(as.Date('2027-03-22'), as.Date('2027-03-26'), 1),  # Spring Break
    # seq(as.Date('2027-05-31'), as.Date('2027-05-31'), 1)   # Memorial Day
  ),
  value    := 1.0]
  
  # add college holidays
  d_college_holidays[date %in% c(
    # --- 2022-23 school year (source: UW-Madison) ---
    seq(as.Date('2023-01-01'), as.Date('2023-01-16'), 1),  # Winter Break
    seq(as.Date('2023-03-13'), as.Date('2023-03-17'), 1),  # Spring Break
    seq(as.Date('2023-05-15'), as.Date('2023-09-05'), 1),  # Summer Break
    
    # --- 2023-24 school year (source: UW-Madison) ---
    seq(as.Date('2023-11-23'), as.Date('2023-11-24'), 1),  # Fall/Thanksgiving Break 
    seq(as.Date('2023-12-22'), as.Date('2024-01-15'), 1),  # Winter Break (extrapolated)
    seq(as.Date('2024-03-25'), as.Date('2024-03-29'), 1),  # Spring Break (extrapolated)
    seq(as.Date('2024-05-13'), as.Date('2024-09-04'), 1),  # Summer Break (extrapolated)
    
    # --- 2024-25 school year (source: UW-Madison) ---
    seq(as.Date('2024-11-28'), as.Date('2024-11-29'), 1),  # Fall/Thanksgiving Break
    seq(as.Date('2024-12-20'), as.Date('2025-01-20'), 1),  # Winter Break
    seq(as.Date('2025-03-24'), as.Date('2025-03-28'), 1),  # Spring Break
    seq(as.Date('2025-06-12'), as.Date('2025-09-02'), 1),  # Summer Break
    
    # --- 2025-26 school year (source: UW-Madison) ---
    seq(as.Date('2025-11-27'), as.Date('2025-11-28'), 1),  # Fall/Thanksgiving Break
    seq(as.Date('2025-12-19'), as.Date('2026-01-19'), 1),  # Winter Break
    seq(as.Date('2026-03-30'), as.Date('2026-04-03'), 1),  # Spring Break
    seq(as.Date('2026-05-11'), as.Date('2026-09-01'), 1),  # Summer Break
    
    # --- 2026-27 school year (source: UW-Madison) ---
    seq(as.Date('2026-11-26'), as.Date('2026-11-27'), 1),  # Fall/Thanksgiving Break
    seq(as.Date('2026-12-18'), as.Date('2026-12-31'), 1)  # Winter Break
  ),
             value    := 1.0]
  
  # School
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
  
  
  # ########################################################### #
  # ##  2a. Contact reductions: school closures              ####
  # ########################################################### #
  # #       * (pre-, primary and secondary school)
  # 
  # # set default school closure
  # # school_dates_non_holiday <- date_all[!date_all %in% dcal_school_closure]
  # # data.table(category = "schools_closed",
  # #            # date     = seq(as.Date('2020-03-14'),as.Date('2020-06-30'),1),
  # #            date     = school_dates_non_holiday,
  # #            value    = 0.0,
  # #            type = 'double',
  # #            age = NA_integer_,
  # #            stringsAsFactors = F
  # # ) -> d_school_closure
  # # 
  # # d_school_closure[date %in% seq(as.Date('2020-03-14'),as.Date('2020-06-30'),1),value := 1.0]
  # # 
  # # tmp_school_closure <- copy(d_school_closure)
  # # tmp_school_closure[,category:='schools_closed']
  # 
  # d_calendar_holiday[date %in% seq(as.Date('2020-03-14'),as.Date('2020-06-30'),1), value := 1.0] 
  # 
  # # set eligible dates for school reopening in May-June 2020
  # d_school_reopening <- seq(as.Date('2020-05-18'),as.Date('2020-06-30'),1)
  # d_school_reopening_wday <- as.POSIXlt(d_school_reopening)$wday
  # 
  # # preschool (reopens 4d/week)
  # d_school_reopening_4d <- d_school_reopening[d_school_reopening_wday %in% 1:4]
  # d_calendar_holiday[date %in% d_school_reopening_4d & age %in% c(0,1,2,6,7),value:=0.5]
  # 
  # # primary school (reopens 2d/week)
  # d_school_reopening_2d <- d_school_reopening[d_school_reopening_wday %in% 4:5]
  # d_school_reopening_2d[1:2] <- d_school_reopening_2d[1:2] - 2 # fix for holidays Thu-Fri in May
  # d_calendar_holiday[date %in% d_school_reopening_2d & age %in% c(11),value:=0.5]
  # 
  # 
  # #secondary school (reopens 1d week)
  # d_school_reopening_1d <- d_school_reopening[d_school_reopening_wday %in% 3]
  # d_calendar_holiday[date %in% d_school_reopening_1d & age %in% c(17),value:=0.5]
  # 
  # # school reopening September 2020
  # # up to primary school
  # d_calendar_holiday[date >= as.Date('2020-09-01') &
  #                      age <= 12 &
  #                    value == 0, value := 0.5]
  # # secondary school
  # d_calendar_holiday[date >= as.Date('2020-05-01') &
  #                      age > 12 & age < 18 &
  #                      value != 1, value := 0.2]
  # # tertiary eduction
  # d_calendar_holiday[date >= as.Date('2020-09-01') &
  #                      age >= 18 &
  #                      value == 0, value := 0.3]
  # 
  # # tertiary eduction: closed from November 1st
  # d_calendar_holiday[date >= as.Date('2020-11-01') &
  #                      age >= 18 , value := 1]

  ########################################### #
  ##  2b. Contact reductions: other        ####
  ########################################### #
  #       * workplace
  #       * community
  #       * household clusters
  
  # # workplace distancing
  # data.table(category = "workplace_distancing",
  #            #date     = seq(as.Date('2020-03-14'),as.Date('2020-05-03'),1),
  #            date     = seq(as.Date(date_start),as.Date(date_end),1),
  #            value    = 0.0,
  #            type     = 'double',
  #            age = NA_integer_,
  #            stringsAsFactors = F
  # ) -> dcal_workplace_distancing
  # dcal_workplace_distancing[date %in% seq(as.Date('2020-03-14'),as.Date('2020-05-03'),1),value := 1.0]
  # 
  # 
  # # community distancing
  # data.table(category = "community_distancing",
  #            #date     = seq(as.Date('2020-03-14'),as.Date('2020-05-24'),1),
  #            date     = seq(as.Date(date_start),as.Date(date_end),1),
  #            value    = 0.0,
  #            type = 'double',
  #            age = NA_integer_,
  #            stringsAsFactors = F
  # ) -> dcal_community_distancing
  # dcal_community_distancing[date %in% seq(as.Date('2020-03-14'),as.Date('2020-05-24'),1),value := 1.0]
  
  # # collectivity distancing
  # data.table(category = "collectivity_distancing",
  #            date     = seq(as.Date(date_start),as.Date(date_end),1),
  #            value    = 0.0,
  #            type = 'double',
  #            age = NA_integer_,
  #            stringsAsFactors = F
  # ) -> dcal_collectivity_distancing
  
  # # household clustering
  # data.table(category = "household_clustering",
  #            date     = seq(as.Date('2020-05-11'),as.Date('2020-08-31'),1),
  #            value    = 1,
  #            type = 'boolean',
  #            age = NA_integer_,
  #            stringsAsFactors = F
  # ) -> dcal_household_clustering
  
  ########################################### #
  ## 3. Imported cases                 ####
  ########################################### #
  
  # data.table(category = "imported_cases",
  #            # date     = seq(as.Date('2020-07-01'),as.Date('2020-08-31'),1),
  #            date     = seq(as.Date(date_start),as.Date(date_end),1),
  #            value    = 0,
  #            type = 'boolean',
  #            age = NA_integer_,
  #            stringsAsFactors = F
  # ) -> dcal_imported_cases
  # dcal_imported_cases[date %in% seq(as.Date('2020-07-01'),as.Date('2020-08-31'),1),value := 1]
  
  ########################################### #
  ##  4. Contact tracing                 ####
  ########################################### #
  
  # data.table(category = "contact_tracing",
  #            date     = seq(as.Date('2020-05-11'),as.Date('2020-08-31'),1),
  #            value    = 1,
  #            type = 'boolean',
  #            age = NA_integer_,
  #            stringsAsFactors = F
  # ) -> dcal_contact_tracing
  
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
  
  plot_calendar(dt_calendar            = d_calendar_all,
                filename_calendar_full = filename_calendar_full,
                show_plots             = show_plots)
  
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

plot_calendar <- function(dt_calendar, filename_calendar_full, show_plots = TRUE){
  
  if(show_plots){
    
    # open pdf stream
    pdf(file=gsub('.csv','.pdf',filename_calendar_full),6,6)
    
    category_opt <- unique(dt_calendar$category)
    par(mfrow=c(3,2))
    
    # check if dt_calendar is data.table
    if(!is.data.table(dt_calendar)){
      dt_calendar <- data.table(dt_calendar)
    }
    
    # make sure that "date" is in date format
    dt_calendar$date <- as.Date(dt_calendar$date)
    
    # x_lim      <- range(dt_calendar$date)
    x_lim      <- c(as.Date('2023-01-01'),max(dt_calendar$date))
    i_cat <- category_opt[2]
    for(i_cat in category_opt){
      plot(x   = dt_calendar[category == i_cat,date],
           y   = dt_calendar[category == i_cat,value],
           xlim = x_lim,
           ylim = range(0,1,dt_calendar$value[dt_calendar$category == i_cat]),
           col  = 1,
           #type='l',
           pch  = 15,
           #lwd=2,
           main = i_cat,
           bty='n',
           xlab = '',
           ylab = unique(dt_calendar[category == i_cat,type]),
           xaxt = 'n'
      )
      add_x_axis(x_lim)
      abline(h=1,lty=3,col='grey')
    }
    
    if("schools_closed" %in% dt_calendar$category){
      
      i_cat <- "schools_closed"

      # plot by age
      plot(x   = dt_calendar[category == i_cat & value == 1,date],
           y   = dt_calendar[category == i_cat & value == 1,age],
           xlim = x_lim,
           ylim = range(0,1,dt_calendar$age,na.rm=T),
           col  = 1,
           pch  = 15,
           main = i_cat,
           bty='n',
           xlab = '',
           ylab = 'age',
           xaxt = 'n'
      )
      add_x_axis(x_lim)
    }
    
    # close pdf stream
    dev.off()
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
adjust_calendar_file <- function(db_category, db_update, file_name, db_age = 'NA', show_plots=FALSE,
                                 bool_singletons = FALSE, erase_category = TRUE){
  
  # file_name fix, exclude '../'
  file_name <- gsub('\\.\\.','\\.',file_name)
  
  # read calendar file
  d_calendar_all <- data.table(read.table(file=file_name,sep=',',header=T))
  
  # make sure the value is stored as double
  d_calendar_all$value <- as.double(d_calendar_all$value)
  
  # adjust age == NA
  d_calendar_all[, age_char := as.character(age)]
  d_calendar_all[is.na(age_char), age_char := 'NA']
  
  d_calendar_categories <- c('general', 
                             'schools_closed', 
                             'workplace_distancing',
                             'collectivity_distancing', 
                             'community_distancing', 
                             'contact_tracing',
                             'contact_survey',
                             'household_clustering', 
                             'imported_cases')
  
  # check category
  if(!db_category %in% d_calendar_categories){
    smd_print("CALENDAR CATEGORY UNKNOWN => STOP CALENDAR ADJUSTMENT")
    smd_print("CALENDAR CATEGORY OPTIONS:", paste0(d_calendar_categories,collapse=', '))
    return(NA)
  }
  
  # remove existing values of this category?
  if(erase_category){
    d_calendar_all <- d_calendar_all[category != db_category,]
  }

  # create data.frame with all unique information to extrapolate
  df_update  <- data.frame(t(db_update))
  df_update  <- unique(df_update)
  
  # define start and end date
  date_out   <- seq(min(as.Date(df_update[,1])),max(as.Date(df_update[,1])),1)
  
  if(length(d_calendar_all$date)>0){
    date_out   <- date_out[date_out<=max(d_calendar_all$date,na.rm = T)]    
  }

  # extrapolate given dates and values
  df_update_full <- list(date = df_update[,1],
                         value = as.numeric(df_update[,2]))
  if(!bool_singletons){
  df_update_full  <- approx(x=as.Date(df_update[,1]),
                            y=df_update[,2],
                            xout = as.Date(date_out),
                            method="linear")
  names(df_update_full) <- c('date','value')
  }
  
  # exclude '0'
  df_update_full$date  <- df_update_full$date[df_update_full$value != 0]
  df_update_full$value <- df_update_full$value[df_update_full$value != 0]
  
  # integrate (new) values in calendar
  for(i_db_age in as.character(db_age)){
    # extrapolate new values
    dcal_new <- data.table(category = db_category,
                           date     = paste(df_update_full$date),
                           value    = df_update_full$value,
                           type     = 'double',
                           age = ifelse(i_db_age == 'NA', NA_integer_,as.numeric(i_db_age)),
                           age_char = i_db_age,
                           stringsAsFactors = F
    ) 
    if(nrow(d_calendar_all)==0){
      d_calendar_all <- dcal_new
    } else {
      # remove existing values for these dates and ages (if any)
      d_calendar_all <- d_calendar_all[!(as.character(date) %in% as.character(df_update_full$date) &
                                           category == db_category &
                                           age_char == i_db_age),]
      
      # add new values
      d_calendar_all <- rbind(d_calendar_all,dcal_new) 
    }
    
  }

  # check
  d_calendar_all[as.character(date) %in% as.character(df_update_full$date) &
                   category == db_category &
                   age_char %in% db_age]
  
  # explore
  plot_calendar(dt_calendar            = d_calendar_all,
                filename_calendar_full = file_name,
                show_plots             = show_plots)
  
  # save as csv 
  write.table(d_calendar_all,
              file = file_name,sep=',',row.names=F,quote=F)
  
}

replace_calendar_value <- function(file_name,db_category,value_orig,value_new,show_plots){
  
  # read calendar file
  d_calendar_all <- data.table(read.table(file=file_name,sep=',',header=T))
  
  # make sure the value is stored as double
  d_calendar_all$value <- as.double(d_calendar_all$value)
  
  d_calendar_all[category == db_category &
                   value == value_orig, value:= value_new]
  
  # explore
  plot_calendar(dt_calendar            = d_calendar_all,
                filename_calendar_full = file_name,
                show_plots             = show_plots)
  
  # save as csv 
  write.table(d_calendar_all,
              file = file_name,sep=',',row.names=F,quote=F)
}

exclude_calendar_category <- function(file_name,db_category,show_plots=FALSE){
  
  # read calendar file
  d_calendar_all <- data.table(read.table(file=file_name,sep=',',header=T))
  
  # remove category (if present)
  d_calendar_all <- d_calendar_all[category != db_category,]
  
  # explore
  plot_calendar(dt_calendar            = d_calendar_all,
                filename_calendar_full = file_name,
                show_plots             = show_plots)
  
  # save as csv 
  write.table(d_calendar_all,
              file = file_name,sep=',',row.names=F,quote=F)
}

# create calendar file comparable to the original lockdown/exit parameter structure
# note: erase_category is a boolean to remove existing values of the given categories 
integrate_parameters_in_calendar <- function(config_exp,
                                             bool_maintain_file_name = FALSE,
                                             erase_category = TRUE){

  # check if social contact survey is activated, but no explicit dates are set
  if(config_exp$event_log_level == "Participants"){
    if(!"contact_survey_dates" %in% names(config_exp) || is.na(config_exp$contact_survey_dates)){
      config_exp$contact_survey_dates <- c_str(as.Date(config_exp$start_date) + 0:(config_exp$num_days-1))
    }
  }
  
  # if there are no time-specific parameters, return original config_exp
  param_calendar <- config_exp[grepl('cnt_reduction_workplace',names(config_exp)) |   # OR colname contains reduction_workplace
                                   grepl('clustering',names(config_exp)) |            # OR colname contains clustering
                                   grepl('imported',names(config_exp)) |              # OR colname contains imported
                                   grepl('survey',names(config_exp)) |
                                   grepl('contact_tracing_date',names(config_exp)) |
                                   grepl('distancing',names(config_exp)) &            # OR colname contains distancing)
                                   !is.na(config_exp)]                                # AND different from NA 
  param_calendar <- unlist(param_calendar)
  param_calendar[is.na(param_calendar)] <- 0

  if(length(param_calendar) == 0 || all(param_calendar == 0)){
    return(config_exp)
  }
  
  # else, modify calendar
  file_name_exp <- config_exp$holidays_file
  
  if(file.exists(file_name_exp)){
    if(!bool_maintain_file_name)
    {
      file_name_new <- smd_file_path(config_exp$output_prefix, gsub('.csv', '_param.csv', basename(config_exp$holidays_file)))
      file.copy(from = file_name_exp,
                to = file_name_new,overwrite = TRUE)
      config_exp$holidays_file <- file_name_new
    }
  } else{
    config_exp$holidays_file <- create_calendar_file(file_name = file_name_exp, show_plots = T)
  }
 
  if('distancing_workplace_ratio' %in% names(config_exp)){
    include_temporal_distancing_factors(db_category    = 'workplace_distancing',
                                        db_values_char = config_exp$distancing_workplace_ratio,
                                        db_delay_char  = config_exp$distancing_workplace_delay,
                                        file_name      = config_exp$holidays_file,
                                        show_plots     = T,
                                        db_dates_char  = config_exp$distancing_workplace_date,
                                        erase_category = erase_category)
  }
  
  if('distancing_school_ratio' %in% names(config_exp)){
    include_temporal_distancing_factors(db_category    = 'schools_closed',
                                        db_values_char = config_exp$distancing_school_ratio,
                                        db_age_char    = c_str(0:25),
                                        db_delay_char  = config_exp$distancing_school_delay,
                                        file_name      = config_exp$holidays_file,
                                        show_plots     = T,
                                        db_dates_char  = config_exp$distancing_school_date,
                                        erase_category = FALSE)
  }
  
  if('distancing_community_ratio' %in% names(config_exp)){
    include_temporal_distancing_factors(db_category    = 'community_distancing',
                                        db_values_char = config_exp$distancing_community_ratio,
                                        db_delay_char  = config_exp$distancing_community_delay,
                                        file_name      = config_exp$holidays_file,
                                        show_plots     = T,
                                        db_dates_char  = config_exp$distancing_community_date,
                                        erase_category = erase_category)
  }
  
  if('distancing_collectivity_ratio' %in% names(config_exp)){
    include_temporal_distancing_factors(db_category    = 'collectivity_distancing',
                                        db_values_char = config_exp$distancing_collectivity_ratio,
                                        db_delay_char  = config_exp$distancing_collectivity_delay,
                                        file_name      = config_exp$holidays_file,
                                        show_plots     = T,
                                        db_dates_char  = config_exp$distancing_collectivity_date,
                                        erase_category = erase_category)
  }
  
  if('imported_cases_number' %in% names(config_exp)){
    include_temporal_distancing_factors(db_category    = 'imported_cases',
                                        db_values_char = config_exp$imported_cases_number,
                                        db_delay_char  = config_exp$imported_cases_delay,
                                        file_name      = config_exp$holidays_file,
                                        show_plots     = T,
                                        db_dates_char  = config_exp$imported_cases_date,
                                        erase_category = erase_category)
  }
  
  if('household_clustering_ratio' %in% names(config_exp)){
    include_temporal_distancing_factors(db_category    = 'household_clustering',
                                        db_values_char = config_exp$household_clustering_ratio,
                                        db_delay_char  = config_exp$household_clustering_delay,
                                        file_name      = config_exp$holidays_file,
                                        show_plots     = T,
                                        db_dates_char  = config_exp$household_clustering_date,
                                        erase_category = erase_category)
  }
  if('contact_survey_dates' %in% names(config_exp)){
    include_temporal_distancing_factors(db_category    = 'contact_survey',
                                        db_values_char = 1,
                                        db_delay_char  = 0,
                                        file_name      = config_exp$holidays_file,
                                        show_plots     = T,
                                        db_dates_char  = config_exp$contact_survey_dates,
                                        bool_singletons= TRUE,
                                        erase_category = erase_category)
  }
  if('contact_tracing_date' %in% names(config_exp)){
    
    include_temporal_distancing_factors(db_category    = 'contact_tracing',
                                        db_values_char = 1,
                                        db_delay_char  = 0,
                                        file_name      = config_exp$holidays_file,
                                        show_plots     = T,
                                        db_dates_char  = config_exp$contact_tracing_date,
                                        bool_singletons= FALSE,
                                        erase_category = erase_category)
  }
  # return list
  return(config_exp)
}

# db_category <- 'distancing_workplace'
# db_values <- seq(0.8,0.9,length=12)
include_temporal_distancing_factors <- function(db_category,db_values_char,
                                                db_age_char = NA,
                                                db_delay_char,
                                                file_name,
                                                show_plots=T,
                                                db_dates_char=NA,
                                                bool_singletons = FALSE,
                                                erase_category  = TRUE){
  
  # check input parameters
  vector_input_param <- c(db_category,db_values_char,db_delay_char,db_dates_char)
  if(!any(is.na(vector_input_param)) & all(vector_input_param != "NA"))
  {
   
  # make sure input parameters are "character" types
  db_values_char <- c_str(db_values_char)
  db_delay_char  <- c_str(db_delay_char)
  db_dates_char  <- paste(db_dates_char)
  
  # split given value string, into numerical values
  db_values <- as.numeric(unlist(strsplit(db_values_char,',')))
  db_delay  <- as.numeric(unlist(strsplit(db_delay_char,',')))
  db_dates  <- as.Date(unlist(strsplit(db_dates_char,',')))
  
  if(!all(is.na(db_age_char))){
    db_age  <- as.numeric(unlist(strsplit(db_age_char,',')))
  } else {
    db_age  <- as.numeric(db_age_char)
  }
  
  if(!bool_singletons){
    # account for delay == 0 by using "db_date-1" and "delay 1" 
    bool_delay_zero <- db_delay == 0
    if(any(bool_delay_zero)){
      db_delay[bool_delay_zero] <- 1
      db_dates[bool_delay_zero] <- db_dates[bool_delay_zero] - 1
    }
   
    # include dates for the delay in compliance
    db_dates  <- c(db_dates[1],db_dates + db_delay,db_dates[-1])
    db_values <- c(0,db_values,db_values[-length(db_values)])
    
    # sort
    db_values <- db_values[order(as.Date(db_dates))]
    db_dates  <- db_dates[order(as.Date(db_dates))]
    
    # add right tail
      db_dates  <- c(db_dates,max(db_dates)+356*3)
      db_values <- c(db_values,db_values[length(db_values)])
  } else {
    if(length(db_values)==1){ db_values <- rep(db_values,length(db_dates))}
  } 
  
  # sort
  db_values <- db_values[order(as.Date(db_dates))]
  db_dates  <- db_dates[order(as.Date(db_dates))]
  
  adjust_calendar_file(db_category = db_category,
                       db_update   = rbind(as.character(db_dates),db_values),
                       db_age      = db_age,
                       file_name   = file_name,
                       show_plots  = TRUE,
                       bool_singletons = bool_singletons,
                       erase_category  = erase_category)
  
  } # end if-clause on is.na
}
