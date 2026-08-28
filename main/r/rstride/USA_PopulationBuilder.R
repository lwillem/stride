############################################################################ #
#
# PREPARE USA-SPECIFIC POPULATION FILE FOR STRIDE
# 
############################################################################ #

#clear workspace
rm(list=ls())

## Load packages
library(haven)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(VGAM)
library(mgcv)
library(arrow)
library(XML)
# library(sf)
# library(sfarrow)
# library(tigris)
library(usmap)

source("~/Documents/Repositories/stride/main/r/rstride/factories/generate-population-community.R")
options(scipen=999)    # turn off scientific notation
# options(scipen=0)      # turn on scientific notation

############################### Synthetic FRED Population ################################

## Manually download data from FRED: https://fred.publichealth.pitt.edu/syn_pops

# MN - Hennepin
# MKE - Milwaukee
# IL - Cook

getFREDdata <- function(state = "WI", county = "Milwaukee", export = FALSE){
  
  print("FRED population files for the desired county must be manually downloaded from 'https://fred.publichealth.pitt.edu/syn_pops'
  and stored in the 'data/population_usa/' folder before running this script.")
  
  if(state == "CT") {       ## Use old data to access old county names for CT
    folder <- fips(state, county = county, data_year = 2021)
    state_name_full <- fips_info(folder, data_year = 2021)$full
    cty_name_full <- fips_info(folder, data_year = 2021)$county
  } else {
    folder <- fips(state, county = county)
    state_name_full <- fips_info(folder)$full
    cty_name_full <- fips_info(folder)$county
  }

  dir <- grep(folder, list.dirs(paste0("data/population_usa/", state)), value = TRUE)
  
  if(length(dir) == 0){
    stop("!!! County data not found !!! \nCheck if population files are stored in 'data/population_usa/' 
         or download county data from 'https://fred.publichealth.pitt.edu/syn_pops'.")
  } else {
    print(paste("Processing population files from:", dir))
    files <- gsub(".txt", "", list.files(dir))
    files <- files[!grepl("METADATA|hospitals", files)]
  }
  
  # import FRED population files
  for(f in files){
    tmp <- read.table(paste0(dir, "/", f, ".txt"), header = TRUE)
    # assign(paste(f, folder, sep = "_"), tmp)
    assign(f, tmp)
  }

  if(min(nchar(households$stcotrbg)) != 12){
    households$stcotrbg <- sprintf("%012.0f", households$stcotrbg)
  }  
  
  # return(list(
  #   gq_people = gq_people,
  #   gq = gq,
  #   households = households,     ## stcotrbg = state, county, tract, blockgroup
  #   people = people,
  #   schools = schools,
  #   workplaces = workplaces,
  #   state = state,
  #   county = county
  # ))
  
  # assign household ids to individuals in people file
  people_data <- inner_join(people, 
                         households, #[wi_025$households$sp_id %in% hh_sample,], 
                         by = c("sp_hh_id" = "sp_id"))
  
  # get census block groups shapefile
  region <- block_groups(state, county, year = 2010) %>% # year aligned with year of FRED population data
    mutate(latitude = as.numeric(INTPTLAT10),
           longitude = as.numeric(INTPTLON10))
  num_bg <- nrow(region)
  district_distance <- matrix(NA,nrow=num_bg,ncol=num_bg)
  i<-1;
  for(i in 1:num_bg){
    dist <- get_distance(region$latitude[i],region$longitude[i],region$latitude,region$longitude)
    district_distance[i,] <- dist
    district_distance[,i] <- dist
  }
  rownames(district_distance) <- region$GEOID
  
  # population community assignment settings
  pop_settings <- data.frame(size                 = nrow(people_data),   # target population size
                             bool_flanders        = FALSE,               # boolean to select flanders
                             bool_teachers        = FALSE,               # boolean to enable "teaching workplaces"
                             bool_census_hh       = FALSE,               # boolean to use the census households, instead of survey data
                             bool_collectivity    = FALSE,               # boolean to generate syntetic nursing homes
                             com_target_size      = 1000,                # target size of the community
                             postfix              = '',                  # to add a tag to the file names
                             # max_age_student      = 23,                # age threshold to participant in school
                             # max_age_teacher      = 60,                # age threshold to be a school teacher
                             rng_seed             = 201909,              # random number generator seed
                             stringsAsFactors     = F) 
  
  # generate population community
  tmp_wp <- generate_population_community(people_data,pop_settings,region,district_distance)
  
  pop_data <- tmp_wp[[1]] %>%
    rename(household_id = sp_hh_id,
           school_id_og = school_id,
           work_id_og = work_id,
           primary_community = com_id,
           secondary_community = com2_id) %>%
    mutate(school_id = dense_rank(school_id_og),
           school_id = ifelse(school_id_og == "X", 0, school_id),
           work_id = dense_rank(work_id_og),
           work_id = ifelse(work_id_og == "X", 0, work_id)
    ) %>%
    select(c("age","household_id","school_id","work_id","primary_community","secondary_community"))
  comm_data <- tmp_wp[[2]]
  
  if(export == TRUE){
    write.table(pop_data,
                paste0("data/pop_US-", state, "-", county, "_c", pop_settings$com_target_size, ".csv"),
                sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)
  }
  
  return(pop_data)
}

# dane <- getFREDdata(state = "WI", county = "Dane")
# mke <- getFREDdata(state = "WI", county = "Milwaukee")
# cook <- getFREDdata(state = "IL", county = "Cook", export = FALSE)
# hennepin <- getFREDdata(state = "MN", county = "Hennepin")

##########################################################################################


