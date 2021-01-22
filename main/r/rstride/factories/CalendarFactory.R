############################################################################ #
#  This file is part of the Stride software. 
#
#  Copyright 2020, Willem L
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
  
  create_calendar_files(21,TRUE)
}


# create calendar files if they do not exist, else re-use them
create_calendar_files <- function(cnt_other_exit_delay = 21, show_plots = FALSE)
{

  filename_calendar_full <- 'data/calendar_belgium_2020_covid19_exit_school_adjusted_universal_import.csv'
  if(file.exists(smd_file_path(filename_calendar_full))){
    return(0) # stop
  }  
  # else, continue
  
########################################### #
## INITIATE DATA                       ####
########################################### #


  
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
           stringsAsFactors = F) -> d_calendar_holiday

summary(d_calendar_holiday)

########################################### #
## 1.b School holidays                 ####
########################################### #

data.table(category = "school_holiday",
           date     = c(seq(as.Date('2019-01-01'),as.Date('2019-01-06'),1), # 2019
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
           value    = 1,
           type = 'boolean',
           age = NA_integer_,
           stringsAsFactors = F
) -> d_school_holidays

data.table(category = "college",
           date     = c(seq(as.Date('2019-01-01'),as.Date('2019-01-06'),1), # 2019
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
           
           value    = 1,
           type = 'boolean',
           age = NA_integer_,
           stringsAsFactors = F
) -> d_college_holidays

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
data.table(category = "schools_closed",
           date     = seq(as.Date('2020-03-14'),as.Date('2020-06-30'),1),
           value    = 1,
           type = 'boolean',
           age = NA_integer_,
           stringsAsFactors = F
) -> d_school_closure

tmp_school_closure <- copy(d_school_closure)
tmp_school_closure[,category:='schools_closed']

dcal_school_closure <- copy(tmp_school_closure[,age:=0])
for(i_age in 1:25){
  dcal_school_closure <- rbind(dcal_school_closure,copy(tmp_school_closure[,age:=i_age]))
}

# set eligible dates for school reopening
d_school_reopening <- seq(as.Date('2020-05-18'),as.Date('2020-06-30'),1)
#d_school_reopening <- d_school_reopening[d_school_reopening != '2020-05-22']
d_school_reopening_wday <- as.POSIXlt(d_school_reopening)$wday

# preschool (reopens 4d/week)
d_school_reopening_4d <- d_school_reopening[d_school_reopening_wday %in% 1:4]
# dcal_preschool_closure <- copy(d_school_closure)
# dcal_preschool_closure[date %in% d_preschool,value:=0]
# dcal_preschool_closure[,category := 'preschool']
dcal_school_closure[date %in% d_school_reopening_4d & age %in% c(0,1,2,6,7),value:=0]

# primary school (reopens 2d/week)
d_school_reopening_2d <- d_school_reopening[d_school_reopening_wday %in% 4:5]
d_school_reopening_2d[1:2] <- d_school_reopening_2d[1:2] - 2 # fix for holidays Thu-Fri in May
# dcal_primary_closure <- copy(d_school_closure)
# dcal_primary_closure[date %in% d_primary,value:=0]
# dcal_primary_closure[,category := 'primary_school']
dcal_school_closure[date %in% d_school_reopening_2d & age %in% c(11),value:=0]


#secondary school (reopens 1d week)
d_school_reopening_1d <- d_school_reopening[d_school_reopening_wday %in% 3]
# d_secondary <- c(d_secondary,d_school_reopening[length(d_school_reopening)]) # add one day in last week
# dcal_secondary_closure <- copy(d_school_closure)
# dcal_secondary_closure[date %in% d_secondary,value:=0]
# dcal_secondary_closure[,category := 'secondary_school']
dcal_school_closure[date %in% d_school_reopening_1d & age %in% c(17),value:=0]


# # college (stay closed)
# dcal_college_closure <- copy(d_school_closure)
# dcal_college_closure[,category := 'college']

########################################### #
##  2b. Contact reductions: other        ####
########################################### #
#       * workplace
#       * community
#       * household clusters

# workplace distancing
data.table(category = "workplace_distancing",
           date     = seq(as.Date('2020-03-14'),as.Date('2020-05-03'),1),
           value    = 1,
           type = 'boolean',
           age = NA_integer_,
           stringsAsFactors = F
) -> dcal_workplace_distancing


# community distancing
data.table(category = "community_distancing",
           date     = seq(as.Date('2020-03-14'),as.Date('2020-05-03')+cnt_other_exit_delay,1),
           value    = 1,
           type = 'boolean',
           age = NA_integer_,
           stringsAsFactors = F
) -> dcal_community_distancing

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

########################################### #
## EXPLORE DATA                         ####
########################################### #

# open pdf stream
#pdf(file='./sim_output/calendar_profile.pdf',6,6)


plot_calendar(d_calendar_all,show_plots)
plot_calendar(d_calendar_holiday,show_plots,b_school_repopening = F)

# close pdf stream
#dev.off()

########################################### #
## SAVE AS CSV	 	         ####
########################################### #

# # format date
# d_calendar_all[,date:=format(date,'%Y-%m-%d')]
# format(d_calendar_all$date,'%Y-%m-%d')

# save as csv (none ==>> dummy)
write.table(d_calendar_holiday[category == 'na',],
            file = smd_file_path('sim_output/holidays_none.csv'),sep=',',row.names=F,quote=F)

# save as csv (default holidays)
write.table(d_calendar_holiday,
            file = 'sim_output/holidays_belgium_2019_2021.csv',sep=',',row.names=F,quote=F)

# save as csv (selection)
write.table(d_calendar_all[!category %in% c('universal_testing'),],
            file = 'sim_output/calendar_belgium_2020_covid19_exit_school_adjusted.csv',sep=',',row.names=F,quote=F)

# save as csv (selection)
write.table(d_calendar_all[!category %in% c('import_cases'),],
            file = 'data/calendar_belgium_2020_covid19_exit_school_adjusted_universal.csv',sep=',',row.names=F,quote=F)

# save as csv (all calendar info)
write.table(d_calendar_all,
            file = filename_calendar_full,sep=',',row.names=F,quote=F)

unique(d_calendar_all$category)

#################################################################### #
## SPECIAL CASES: COVID19 school reopening May-June 2020	 	      ####
#################################################################### #

plot_calendar(d_calendar_all,show_plots)

# make copy
d_calendar_exit_school <- copy(d_calendar_all)

# undo all school reopenings
d_calendar_exit_school[value==0,value := 1]

# no school reopening
plot_calendar(d_calendar_exit_school,show_plots)
write.table(d_calendar_exit_school,
            file = 'data/calendar_belgium_2020_covid19_exit_schoolsclosed.csv',sep=',',row.names=F,quote=F)

# preschool
d_calendar_exit_subset <- copy(d_calendar_exit_school)
d_calendar_exit_subset[date %in% d_school_reopening & age %in% 0:5, value := 0]
write.table(d_calendar_exit_subset,
            file = 'data/calendar_belgium_2020_covid19_may_preschool.csv',sep=',',row.names=F,quote=F)
plot_calendar(d_calendar_exit_subset,show_plots)


# primary school
d_calendar_exit_subset <- copy(d_calendar_exit_school)
d_calendar_exit_subset[date %in% d_school_reopening & age %in% 6:11, value := 0]
write.table(d_calendar_exit_subset,
            file = 'data/calendar_belgium_2020_covid19_may_primary_school.csv',sep=',',row.names=F,quote=F)
plot_calendar(d_calendar_exit_subset,show_plots)

# secondary school
d_calendar_exit_subset <- copy(d_calendar_exit_school)
d_calendar_exit_subset[date %in% d_school_reopening & age %in% 12:17, value := 0]
write.table(d_calendar_exit_subset,
            file = 'data/calendar_belgium_2020_covid19_may_secondary_school.csv',sep=',',row.names=F,quote=F)
plot_calendar(d_calendar_exit_subset,show_plots)

# K6 schools
d_calendar_exit_subset <- copy(d_calendar_exit_school)
d_calendar_exit_subset[date %in% d_school_reopening & age %in% 0:11, value := 0]
write.table(d_calendar_exit_subset,
            file = 'data/calendar_belgium_2020_covid19_may_k6school.csv',sep=',',row.names=F,quote=F)
plot_calendar(d_calendar_exit_subset,show_plots)

#K12 school
d_calendar_exit_subset <- copy(d_calendar_exit_school)
d_calendar_exit_subset[date %in% d_school_reopening & age %in% 0:17, value := 0]
write.table(d_calendar_exit_subset,
            file = 'data/calendar_belgium_2020_covid19_may_k12school.csv',sep=',',row.names=F,quote=F)
plot_calendar(d_calendar_exit_subset,show_plots)

# pre-, primary and secondary fully open on selection of days
# note: based on the age-specific program, but by school category
d_calendar_exit_subset <- copy(d_calendar_exit_school)
d_calendar_exit_subset[date %in% d_school_reopening_4d & age %in% 0:5, value := 0]
d_calendar_exit_subset[date %in% d_school_reopening_2d & age %in% 6:11, value := 0]
d_calendar_exit_subset[date %in% d_school_reopening_1d & age %in% 12:17, value := 0]
write.table(d_calendar_exit_subset,
            file = 'data/calendar_belgium_2020_covid19_exit_schoolcategory_adjusted.csv',sep=',',row.names=F,quote=F)
plot_calendar(d_calendar_exit_subset,show_plots)

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
    i_cat <- category_opt[1]
    for(i_cat in category_opt){
      plot(x   = dt_calendar[category == i_cat,date],
           y   = dt_calendar[category == i_cat,value],
           xlim = x_lim,
           ylim = range(0,1,dt_calendar$value),
           col  = 1,
           pch  = 15,
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
      add_x_axis(x_lim)
    }
  }
}

