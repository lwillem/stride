############################################################################ #
#  This file is part of the Stride software. 
#
#  Copyright 2020, Willem L, Libin P
############################################################################ #
#
# TO CREATE CALENDAR FILE(S) FOR UNIVERSAL TESTING SIMULATIONS
#
# CONTAINING:
# 1. Public and school holidays
# 2. Contact reductions (covid19)
#       * pre-, primairy and secondary school
#       * workplace
#       * community
#       * household clusters     ==>> NOT USED HERE
# 3. Imported cases (covid19)
# 4. Contact tracing (covid19)   ==>> NOT USED HERE
# 5. Universal testing (covid19)
#
# GOAL(S): 
# - to switch from lockdown to universal testing
# - what happens if cases are imported when local transmission is controlled  
#
# NOTES: 
# - In this script, all variables with a name starting with "dcal_" are 
#   automatically merged.
# - Creates one csv calendar file with and one without imported cases
#
############################################################################ #
## code for debugging ----
if(0==1){
  
  # specify delay to import cases
  delay_import_cases <- 8*7
  
  # inspect calendar files
  show_plots <- TRUE
  
  # create calender with switch on the 1th of May
  date_policy_switch <- "2020-05-01"
  create_calenders_universal_testing(date_policy_switch,delay_import_cases)

  # create calender with switch on the 1th of July
  create_calenders_universal_testing("2020-07-01",delay_import_cases)
}

create_calenders_universal_testing <- function(date_policy_switch, 
                                               delay_import_cases,
                                               show_plots = FALSE)
{
  
  ########################################### #
  ## INITIATE DATA                       ####
  ########################################### #
  
  ## default value in C++ CALENDAR vectors ==>> 0
  
  ## make sure the given date to switch policy is of type "date"
  date_policy_switch <- as.Date(date_policy_switch)
  
  ## set lockdown period
  dates_lockdown <- seq(as.Date('2020-03-14'),date_policy_switch-1,1)
  
  ## set lockdown period
  dates_universal_testing <- seq(date_policy_switch,as.Date('2021-12-31'),1) # use a dummy end (has no impact on the simulation)
  
  ## set dates for imported cases
  dates_import_cases <- dates_universal_testing + delay_import_cases

  ########################################### #
  ## 1.a Public holidays                 ####
  ########################################### #
  
  data.table(category = "general",
             date     = as.Date(c('2020-01-01','2020-04-13','2020-05-01','2020-05-21','2020-06-01', # 2020
                                  '2020-07-21','2020-08-15','2020-11-01','2020-11-11','2020-12-25')),
             value    = 1,
             type = 'boolean',
             age = NA_integer_,
             stringsAsFactors = F) -> d_calendar_holiday
  
  summary(d_calendar_holiday)
  
  ########################################### #
  ## 1.b School holidays                 ####
  ########################################### #
  
  data.table(category = "school_holiday",
             date     = c(seq(as.Date('2020-01-01'),as.Date('2020-01-05'),1), # 2020
                          seq(as.Date('2020-02-24'),as.Date('2020-02-29'),1),
                          seq(as.Date('2020-04-06'),as.Date('2020-04-19'),1),
                          seq(as.Date('2020-07-01'),as.Date('2020-08-31'),1),
                          seq(as.Date('2020-11-02'),as.Date('2020-11-08'),1),
                          seq(as.Date('2020-12-21'),as.Date('2020-12-31'),1)),
             value    = 1,
             type = 'boolean',
             age = NA_integer_,
             stringsAsFactors = F
  ) -> d_school_holidays
  
  data.table(category = "college",
             date     = c(seq(as.Date('2020-01-01'),as.Date('2020-01-05'),1), # 2020
                          seq(as.Date('2020-02-24'),as.Date('2020-02-29'),1),
                          seq(as.Date('2020-04-06'),as.Date('2020-04-19'),1),
                          seq(as.Date('2020-07-01'),as.Date('2020-09-20'),1),# summer break untill September, 20
                          #seq(as.Date('2020-11-02'),as.Date('2020-11-08'),1), # no fall break
                          seq(as.Date('2020-12-21'),as.Date('2020-12-31'),1)),
             value    = 1,
             type = 'boolean',
             age = NA_integer_,
             stringsAsFactors = F
  ) -> d_college_holidays
  
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
  
  ########################################################### #
  ##  2a. Contact reductions: school closures              ####
  ########################################################### #
  #       * (pre-, primary and secondary school)
  
  # set default school closure
  sel_ages <- 0:17
  data.table(category = "schools_closed",
             date     = rep(dates_lockdown,length(sel_ages)),
             value    = 1,
             type     = 'boolean',
             age      = rep(sel_ages,each=length(dates_lockdown)),
             stringsAsFactors = F
  ) -> dcal_school_closure
  
  #secondary school
  dcal_secondary_closure <- copy(dcal_school_closure)
  dcal_secondary_closure[,category := 'secondary_school']
  
  # college (default remains closed)
  sel_ages  <- 18:25
  sel_dates <- c(dates_lockdown,dates_universal_testing)
  data.table(category = "schools_closed",
             date     = rep(sel_dates,length(sel_ages)),
             value    = 1,
             type = 'boolean',
             age = rep(sel_ages,each=length(sel_dates)),
             stringsAsFactors = F
  ) -> dcal_college_closure

  
  ########################################### #
  ##  2b. Contact reductions: other        ####
  ########################################### #
  #       * workplace
  #       * community
  #       * household clusters (not included here)
  
  # workplace distancing
  data.table(category = "workplace_distancing",
             date     = dates_lockdown,
             value    = 1,
             type = 'boolean',
             age = NA_integer_,
             stringsAsFactors = F
  ) -> dcal_workplace_distancing
  
  
  # workplace distancing
  data.table(category = "community_distancing",
             date     = dates_lockdown,
             value    = 1,
             type = 'boolean',
             age = NA_integer_,
             stringsAsFactors = F
  ) -> dcal_community_distancing
  
  
  ########################################### #
  ## 3. Imported cases                 ####
  ########################################### #
  
  data.table(category = "imported_cases",
               date     = dates_import_cases,
               value    = 1,
               type = 'boolean',
               age = NA_integer_,
               stringsAsFactors = F
  ) -> dcal_imported_cases
  
  
  ########################################### #
  ##  4. Contact tracing                 ####
  ########################################### #
  
  ########################################### #
  ## 5. Universal testing                 ####
  ########################################### #
  
  data.table(category = "universal_testing",
             date     = dates_universal_testing,
             value    = 1,
             type = 'boolean',
             age = NA_integer_,
             stringsAsFactors = F
  ) -> dcal_universal_testing
  
  
  ########################################### #
  ## MERGE HOLIDAYS & OTHER CALENDAR ITEMS ####
  ########################################### #
  
  ## make sure, no holidays are present during "universal testing"
  d_calendar_holiday <- d_calendar_holiday[!d_calendar_holiday$date %in% dates_universal_testing,]
  
  # get 'dcal_*' variables
  opt_other <- ls(pattern='dcal_')
  
  # combine all 'dcal_*' variable
  d_calendar_all <- foreach(i_other  = opt_other,
                      .init    = d_calendar_holiday,
                      .combine = 'rbind') %do% {
                        get(i_other)
                      } 
  unique(d_calendar_all$category)
  
  ########################################### #
  ## EXPLORE DATA                         ####
  ########################################### #
  
  plot_calendar(d_calendar_all,show_plots)

  ########################################### #
  ## SAVE AS XML AND CSV	 	         ####
  ########################################### #
  
  # get "day of the year" the policy switch takes place
  day_index_policy_switch  <- strftime(date_policy_switch, format = "%j")
  
  # get "day of the year" the seeding of infected cases (re)starts
  day_index_import_cases   <- strftime(min(dates_import_cases), format = "%j")
  
  # create standardised names
  filename_calendar        <- paste0('calendar_belgium_covid19_universaltest_d',day_index_policy_switch,'.csv')
  filename_calendar_import <- gsub('.csv',paste0('_import_d',day_index_import_cases,'.csv'),filename_calendar)
  filename_calendar_import_college <- gsub('.csv',paste0('_import_college_d',day_index_import_cases,'.csv'),filename_calendar)
  
  # save all calendar info as csv
  write.table(d_calendar_all,
              file = smd_file_path('data',filename_calendar_import),
              sep=',',row.names=F,quote=F)
  
  # save subset (no imported cases) as csv
  write.table(d_calendar_all[!category %in% c('import_cases'),],
              file = smd_file_path('data',filename_calendar),
              sep=',',row.names=F,quote=F)

  # save subset (college reopens)
  sel_ages <- unique(dcal_college_closure[date %in% dates_universal_testing,c(age)])
  d_calendar_all[date %in% dates_universal_testing & category == "schools_closed" & age %in% sel_ages,value:=0]
  write.table(d_calendar_all,
              file = smd_file_path('data',filename_calendar_import_college),
              sep=',',row.names=F,quote=F)
  plot_calendar(d_calendar_all)
  
  smd_print("CREATED", smd_file_path('data',filename_calendar))
  smd_print("CREATED", smd_file_path('data',filename_calendar_import))
  smd_print("CREATED", smd_file_path('data',filename_calendar_import_college))
  
  # return filenames
  return(c(filename_calendar_import,filename_calendar,filename_calendar_import_college))
  
} # end function



