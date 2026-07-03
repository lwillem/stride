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

source("bin/rstride/factories/generate-population-community.R")
options(scipen=999)    # turn off scientific notation
# options(scipen=0)      # turn on scientific notation

############################### Synthetic FRED Population ################################

## Manually download data from FRED: https://fred.publichealth.pitt.edu/syn_pops

# MN - Hennepin
# MKE - Milwaukee
# IL - Cook

getFREDdata <- function(state = "WI", county = "Milwaukee"){
  
  folder <- fips(state, county = county)
  state_name_full <- fips_info(folder)$full
  cty_name_full <- fips_info(folder)$county
  
  dir <- grep(folder, list.dirs(paste0("data/population_usa/", state)), value = TRUE)
  
  files <- gsub(".txt", "", list.files(dir))
  files <- files[!grepl("METADATA|hospitals", files)]
  
  for(f in files){
    tmp <- read.table(paste0(dir, "/", f, ".txt"), header = TRUE)
    # assign(paste(f, folder, sep = "_"), tmp)
    assign(f, tmp)
  }
  
  return(list(
    gq_people = gq_people,
    gq = gq,
    households = households,     ## stcotrbg = state, county, tract, blockgroup
    people = people,
    schools = schools,
    workplaces = workplaces,
    state = state,
    county = county
  ))
}

dane <- getFREDdata(state = "WI", county = "Dane")
mke <- getFREDdata(state = "WI", county = "Milwaukee")
cook <- getFREDdata(state = "IL", county = "Cook")
hennepin <- getFREDdata(state = "MN", county = "Hennepin")

##########################################################################################

pop_data <- inner_join(mke$people, 
                       mke$households, #[wi_025$households$sp_id %in% hh_sample,], 
                       by = c("sp_hh_id" = "sp_id"))

## By Block Group

region <- block_groups(state = mke$state, county = mke$county, year = 2010) %>%
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

# SETTINGS
pop_settings <- data.frame(size                 = nrow(pop_data),      # target population size
                           bool_flanders        = FALSE,               # boolean to select flanders
                           bool_teachers        = FALSE,               # boolean to enable "teaching workplaces"
                           bool_census_hh       = FALSE,               # boolean to use the census households, instead of survey data
                           bool_collectivity    = FALSE,               # boolean to generate syntetic nursing homes
                           com_target_size      = 1000,                 # target size of the community
                           postfix              = '',                  # to add a tag to the file names
                           # max_age_student      = 23,                  # age threshold to participant in school
                           # max_age_teacher      = 60,                  # age threshold to be a school teacher
                           rng_seed             = 201909,              # random number generator seed
                           stringsAsFactors     = F)   

tmp_wp <- generate_population_community(pop_data,pop_settings,region)
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
              "data/pop_US-WI-MKE_c1000.csv",
              sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)
}

