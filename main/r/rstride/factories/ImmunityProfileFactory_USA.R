library(usmap)

# # load data
# us_imm <- read.table("~/Documents/Repositories/measles/data/immunity/Kiang-et-al_measles-immunity-2024.txt", 
#                      header = TRUE)
# 
# rownames(us_imm) <- us_imm$State
# us_imm$State <- NULL
# 
# ## Include Child VaxView immunity estimates from newborns born in 2022 in Kiang immunity estimates
# child_vax <- read.table("~/Documents/Repositories/measles/data/immunity/Vaccination_Coverage_byState_year2022.txt", 
#                         header = TRUE,
#                         sep = "\t") %>%
#   inner_join(fips_info(), by = c("Geography" = "full")) %>%
#   dplyr::select(c(abbr, `Estimate....`, Dimension)) %>%
#   group_by(abbr) %>%
#   spread(value = `Estimate....`, key = Dimension)
# 
# immunity_comb <- inner_join(child_vax, us_imm, by = c("abbr"= "State"))
# ## Assume 13months vax coverage for children <1 year old
# colnames(immunity_comb) <- c("State", "13_months", "1_year", "2_year", "3_years", "4_years",
#                              "5_9_years", "10_14_years", "15_19_years", "20_24_years", "25_29_years",
#                              "30_34_years", "35_39_years", "40_44_years", "45_49_years", "50_54_years",
#                              "55_59_years", "60_64_years", "65_69_years", "70_74_years", "75_79_years", 
#                              "80_85_years", "86_years")
# 
# write.table(immunity_comb, "~/Documents/Repositories/measles/data/immunity/Kiang_ChildVaxView_Immunity_Combined_clean.txt", 
#             sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

get_immunity_usa <- function(state_abbr, max_age = 100){
  # load data
  imm_comb <- read.table("~/Documents/Repositories/measles/data/immunity/Kiang_ChildVaxView_Immunity_Combined_clean.txt", 
                       header = TRUE)
  rownames(imm_comb) <- imm_comb$State
  imm_comb$State <- NULL
  
  # get state data 
  if(nchar(state_abbr) > 2){
    stop("!!! State Abbreviation Not Found !!!")
  } else {
    state_imm <- imm_comb[state_abbr, ]
  }
  
  # expand age groups to max age
  state_imm_exp <- c()
  for(nm in names(state_imm)) {
    nums <- as.integer(unlist(regmatches(nm, gregexpr("\\d+", nm))))
    ages <- if(length(nums) == 1) nums[1] else nums[1]:nums[2]  
    state_imm_exp <- c(state_imm_exp, rep(state_imm[[nm]], length(ages)))
  }
  
  num_age <- max_age+1
  state_imm_exp <- c(state_imm_exp,rep(last(state_imm_exp),num_age-length(state_imm_exp)))
  
  state_imm_exp

  # convert to ratio
  immunity_profile <- state_imm_exp/100
  
  # adjust infant immunity
  # children 6mons+ are eligible for 1 dose vaccine; 
  # children < 6mos are expected to have immunity from mother
  # immunity_profile[1] <- 1/2
  
  # get suscetibility =  1 - immunity
  susceptiblilty_profile <- 1-immunity_profile
  
  # explore
  plot(susceptiblilty_profile,ylim=0:1,type='l',lwd=7,ylab='susceptibility',xlab='age')
  plot(immunity_profile,ylim=0:1,type='l',lwd=7,ylab='immunity',xlab='age')
  
  ############################################
  ## SAVE AS XML  	 	                      ##
  ############################################
  
  # add age group as column names
  names(immunity_profile) <- paste0('age',0:max_age)
  
  # add info on data source and manipulation
  immunity_data <- unlist(list(data_source = 'Kiang_ChildVaxView_Immunity_Combined_clean.txt',
                               data_manipulation = "mean by age",
                               round(immunity_profile,digits=4)))
  
  # save as xml
  .rstride$save_config_xml(immunity_data,'immunity',paste0('immunity_measles_', state_abbr, ".xml"))
}
