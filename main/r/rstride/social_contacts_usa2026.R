############################################################################ #
#  This file is part of the Stride Population software. 
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
#  Copyright 2026, Willem L
############################################################################ #
#
# PREPARE USA-SPECIFIC SOCIAL CONTACT DATA FOR STRIDE
# 
############################################################################ #

#clear workspace
rm(list=ls())

## set working director (or open RStudio with this script)
# setwd("C:/User/path/to/the/r-project/folder") ## WINDOWS
# setwd("/Users/path/to/the/r-project/folder")        ## MAC
# library(here); setwd(here())                     ## set the wd at the location of the 'project-file'

# load 'contactdata' package
suppressPackageStartupMessages(library('contactdata'))

# # load simid.rtools package
# require(devtools)
# devtools::install_github("lwillem/simid_rtools",force=F,quiet=T)
# #devtools::uninstall(simid.rtools)
# library('simid.rtools')
# 
# # update description file, if required
# smd_update_description_file()

# load FRED population builder
source("~/Documents/Repositories/stride/main/r/rstride/USA_PopulationBuilder.R")

# select country with ISO2 code
sel_country <- 'US'

# set state and county
state <- "WI"
county <- "Dane"
export <- FALSE

# select population file
# pop_file <- 'data/pop_US-WI-MKE_c1000.csv'

# create run tag using the current time if use_date_prefix == TRUE
run_tag <- format(Sys.time(), format="%Y%m%d_%H%M%S")

# add dir_postfix
# run_tag <- paste0(run_tag,dir_postfix)

# set output tag
cdata_tag    <- paste0('contact_matrix_usa_conditional')

# set output directory
output_dir <- smd_file_path('sim_output')
folder_name <- file.path(output_dir, paste0(run_tag, "_", state, "-", county, "_conditional_social_contacts"))
if(!dir.exists(folder_name) & export == TRUE) {
  dir.create(folder_name, recursive = TRUE)
}
file_name  <- file.path(folder_name,cdata_tag)

# define help function to load social data from Prem et al. by location
get_cnt_data <- function(location, country) {
   cnt_matrix <- contact_matrix(country = country,
                                   location = location,
                                   #geographic_setting = c("all"),
                                   data_source = c("2020"))
    cnt_count <- rowSums(cnt_matrix)
    all_ages <- 0:94        # convert matrix to vector                     
    approx(x = seq(0,75,5), # expand matrix age groups using linear interpolation
           y = cnt_count,
           xout = all_ages,
           method = 'constant',
           rule = 2:2)$y
}

############################## #
# LOAD AND PROCESS DATA ####
############################## #

# get data by location ----
cnt_all    <- get_cnt_data("all", country = sel_country)
cnt_home   <- get_cnt_data("home", country = sel_country)
cnt_workplace   <- get_cnt_data("work", country = sel_country)#*1.2
cnt_sm_workplace   <- get_cnt_data("work", country = sel_country)*1.2
cnt_school <- get_cnt_data("school", country = sel_country)
cnt_other  <- get_cnt_data("other", country = sel_country)

# check
cnt_all_list <- contact_df_countries(countries = sel_country,
                                     location = "all",
                                     geographic_setting = c("all"),
                                     data_source = c("2020")
                )
sum(cnt_all_list$contact[cnt_all_list$age_from == "75_80"])
cnt_all[76] # note: index 1 is age age 0

# additional contacts ----
# e.g. if contact data assume school contacts for age i, but population data assumes no school enrolment
# e.g. if the number of contacts at home outnumber the number of household members

# define vector to accumulate "additional contact"
cnt_additional <- cnt_all * 0

# demography data ----
# load file
pop_usa <- getFREDdata(state = state, county = county, export = FALSE) #read.csv(pop_file,header = T, sep = ',')

if(export == TRUE){
  write.table(pop_usa, file.path(folder_name, paste0(run_tag, "_population_", state, "-", county, ".csv")),
              sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)
}

# rename work_id to workplace_id
names(pop_usa) <- gsub('work','workplace',names(pop_usa))

# set 'missing' to NA
pop_usa$school_id[pop_usa$school_id == 0] <- NA
pop_usa$workplace_id[pop_usa$workplace_id == 0] <- NA

# define age breaks for the age distribution
breaks_ages <- (0:95) -0.5

# age distribution
age_counts <- hist(pop_usa$age, breaks = breaks_ages, plot = FALSE)$counts

# school enrolment ----
age_counts_school <- hist(pop_usa$age[!is.na(pop_usa$school_id)], breaks = breaks_ages, plot = FALSE)$counts
age_distr_school <- age_counts_school / age_counts
age_distr_school[is.na(age_distr_school)] <- 0

# calculate number of contacts conditional on being at school
cnt_school_conditional      <- cnt_school * (7/5) # school is open 5 out of 7
cnt_school_conditional[3:5] <- cnt_school_conditional[6] # extrapolate behaviour of age 5 to age 2 to 4

# set contact rates for ages not enrolled in school to zero
cnt_school_conditional[age_distr_school == 0] <- 0

# define school contacts for ages not at school as "additional"
cnt_additional     <- cnt_school * (age_distr_school == 0)
cnt_additional[20] <- cnt_additional[21] # adjust artefact for age 19 (not enrolled in US data, but high number of contacts observed)

# employment ----
age_counts_workplace <- hist(pop_usa$age[!is.na(pop_usa$workplace_id)], breaks = breaks_ages, plot = FALSE)$counts
age_distr_workplace <- age_counts_workplace / age_counts
age_distr_workplace[is.na(age_distr_workplace)] <- 0

################################################
######## Account for small workplaces

tmp <- pop_usa %>% group_by(workplace_id) %>% summarize(n= n()) %>% subset(n <= 7)
age_counts_sm_workplace <- hist(pop_usa$age[pop_usa$workplace_id %in% tmp$workplace_id], breaks = breaks_ages, plot = FALSE)$counts

age_distr_workplace_avg <- (age_counts_workplace-age_counts_sm_workplace) / age_counts
age_distr_workplace_avg[is.na(age_distr_workplace_avg)] <- 0
age_distr_workplace_sm <- age_counts_sm_workplace / age_counts
age_distr_workplace_sm[is.na(age_distr_workplace_sm)] <- 0

workplace_ages <- 18:69
workplace_ages_select <- 30:49
cnt_sm_workplace_conditional <- cnt_sm_workplace * 0 # start with zero's
cnt_sm_workplace_conditional[workplace_ages + 1]  <- mean(cnt_sm_workplace[workplace_ages_select + 1]) # index = age + 1
# calculate number of contacts conditional on being at work
cnt_sm_workplace_conditional <- cnt_sm_workplace_conditional * (7/5)  # account for working 5 days out of 7
cnt_sm_workplace_conditional <- cnt_sm_workplace_conditional  / mean(age_distr_workplace[workplace_ages_select+1]) # account for employment rate

################################################################################

# define conditional number of contacts for all ages as the average of a selection of the the (most) active population
workplace_ages <- 18:69
workplace_ages_select <- 30:49
cnt_workplace_conditional <- cnt_workplace * 0 # start with zero's
cnt_workplace_conditional[workplace_ages + 1]  <- mean(cnt_workplace[workplace_ages_select + 1]) # index = age + 1

# calculate number of contacts conditional on being at work
cnt_workplace_conditional <- cnt_workplace_conditional * (7/5)  # account for working 5 days out of 7

cnt_workplace_conditional <- cnt_workplace_conditional / mean(age_distr_workplace[workplace_ages_select+1]) # account for employment rate

# optional: increase rates to account for small workplaces (n < mean number of contacts)?
workplace_size_count <- table(table(pop_usa$workplace_id))
# number of people in workplaces with ≤7 people
num_people_workplace_leq7 <- sum(workplace_size_count[1:7] *  1:7)
# proportional to number of workers
num_people_workplace_leq7 / sum(!is.na(pop_usa$workplace_id)) 

# household contacts ----
# define household sizes
hh_size <- data.frame(table(pop_usa$household_id))
names(hh_size) <- c('household_id','household_size')

# add household info to population matrix
pop_usa_edit <- merge(pop_usa,hh_size)

# count number of households by member age and size
num_hh_age_size <- table(pop_usa_edit$age,pop_usa_edit$household_size) # CHECK THIS
num_age <- rowSums(num_hh_age_size)

# define matrix to represent nubmer of hh contacts by hh size
mat_cnt_size <-  matrix(rep(1:ncol(num_hh_age_size) - 1, length(num_age)), 
                       ncol = ncol(num_hh_age_size), 
                       byrow = T)

# assume fully connected households, get total number of contacts by age
hh_cnt_size <- num_hh_age_size * mat_cnt_size

# get mean number of contacts at home by age if household is fully connected
cnt_home_fully_connected <- rowSums(hh_cnt_size) / num_age

# impute/approximate missing ages
cnt_home_fully_connected <- approx(names(cnt_home_fully_connected),
                                   cnt_home_fully_connected,
                                   0:(length(cnt_home)-1), 
                                   rule = 2)$y  ## For ages > the age group available in population data, copy value from previous age

# set contacts with non-household-members as "additional"
cnt_additional <- cnt_additional +  (cnt_home - cnt_home_fully_connected)

# community
# start from cnt_other and add "additional"
cnt_other_adj <- cnt_other + cnt_additional

# explore ----
# define function to explore (un)conditional contact rates
plot_conditional_contacts <- function(cnt_orig, cnt_conditional, pop_fraction, plot_main, 
                                      state = "", county = "", xlim = c(0,95)){
  
  # define y_limit
  ylim <- c(0, max(c(cnt_orig,cnt_conditional)) * (4/3))
  
  # set margin
  par(mar=c(5,5,2,5))
  
  # plot (un)conditional contact rates
  plot(cnt_conditional,xlim=xlim, 
       ylim = ylim,
       main = paste0(plot_main, " - ", state, ", ", county), 
       xlab = "age", 
       ylab="mean number of contacts")
  points(cnt_orig,col=2)
  legend('topleft',c('unconditional', 'conditional'), fill = 2:1)
  
  # set population fraction scaling factor
  if(!any(is.na(pop_fraction))){
    frac_schale <- max(cnt_conditional) * 3/4
    lines(pop_fraction*frac_schale, col = 4)
    abline(h= frac_schale, lty = 3, col = 4)
    axis(4,seq(0,frac_schale,length.out = 11 ) ,labels = seq(0,10,1)/10, las = 2, col = 4)
    mtext(paste('fraction enrolled in', plot_main), side = 4, padj = 4, col = 4, adj = 0.1)
  }
}

# explore school contacts: conditional and unconditional

# get file name with path
file_name_path <- file.path(folder_name, paste0(run_tag, "_social_contacts_", state, "-", county, "_plots.pdf"))

# check extension and add if not present
if(!grepl('.pdf',file_name_path)){
  file_name_path <- paste0(file_name_path,'.pdf')
}

# open pdf stream
pdf(file_name_path)

plot_conditional_contacts(cnt_school, cnt_school_conditional, age_distr_school, 'school', state = state, county = county, xlim = c(0,22))
plot_conditional_contacts(cnt_workplace, cnt_workplace_conditional, age_distr_workplace, 'workplace', state = state, county = county)
plot_conditional_contacts(cnt_home, cnt_home_fully_connected, NA, 'household', state = state, county = county)
plot_conditional_contacts(cnt_other, cnt_other_adj, NA, 'other', state = state, county = county)
plot_conditional_contacts(cnt_all, cnt_all, NA, 'total', state = state, county = county)

# close pdf stream
dev.off()

###############################
# STORE AS LIST FOR R ####
###############################

# start with info on data and methods
cnt_data_meta <- list(data_source = 'USA social contact based on Prem et al (2017)',
                        method      = "Reported average number of contacts by age group, conditional on presence",
                        author      = Sys.info()['user'],
                        date        = format(Sys.time())
                        )

# add social contact data
social_cnt_data                     <- cnt_data_meta
social_cnt_data$regular_weekday     <- cnt_all
social_cnt_data$regular_weekend     <- cnt_all
social_cnt_data$household           <- cnt_home_fully_connected
social_cnt_data$school              <- cnt_school_conditional
social_cnt_data$workplace           <- cnt_workplace_conditional
social_cnt_data$community_weekday   <- cnt_other_adj
social_cnt_data$community_weekend   <- cnt_other_adj

save(social_cnt_data,file=paste0(file_name,'.RData'))

###############################
## STORE AS XML FOR STRIDE   ##
###############################
library('XML')

# create xml prefix
xml_prefix <- paste0(' This file is part of the Stride software [', format(Sys.time()), ']')


cnt_matrices_lib <- list(regular_weekday     = cnt_all,
                         regular_weekend     = cnt_all,
                         household           = cnt_home_fully_connected,
                         school              = cnt_school_conditional,
                         workplace           = cnt_workplace_conditional,
                         community_weekday   = cnt_other_adj,
                         community_weekend   = cnt_other_adj
                         )
# setup XML doc (to add prefix)
xml_doc = newXMLDoc()

cnt_matrix_xml  <- newXMLNode("matrices", doc = xml_doc)
cnt_matrix_meta <- newXMLNode("metadata", parent = cnt_matrix_xml)
smd_listToXML(cnt_matrix_meta,cnt_data_meta)

i_context <- 1
for(i_context in 1:length(cnt_matrices_lib))
{

  # extract data
  cnt_context_name   <- names(cnt_matrices_lib)[i_context]
  cnt_context_values <- cnt_matrices_lib[[i_context]]

  print(cnt_context_name)

  # add data to XML
  cnt_context <- newXMLNode(cnt_context_name, parent=cnt_matrix_xml)
  for(i in 1:length(cnt_context_values)){

    participant <- newXMLNode("participant",parent=cnt_context)

    part_age  <- newXMLNode("age",parent=participant)
    xmlValue(part_age) <- paste(i)

    contacts  <- newXMLNode("contacts",parent=participant)

    # for(j in 1:ncol(survey_mij)){
       contact        <- newXMLNode("contact",parent=contacts)
       age            <- newXMLNode("age", parent=contact);
       xmlValue(age)  <- 'all'
       rate           <- newXMLNode("rate", parent=contact);
       xmlValue(rate) <- paste(cnt_context_values[i])
    # }
  }
}

# create filename for xml output
out_filename <- paste0(file_name,'.xml')

# xml prefix
xml_prefix <- paste0(' This file is part of the Stride software [', format(Sys.time()), ']')

# save as XML,
# note: if we use an XMLdoc to include prefix, the line break dissapears...
# fix: http://r.789695.n4.nabble.com/saveXML-prefix-argument-td4678407.html
cat(saveXML(xml_doc, indent = TRUE, prefix = newXMLCommentNode(xml_prefix)),  file = out_filename)
print(out_filename)

