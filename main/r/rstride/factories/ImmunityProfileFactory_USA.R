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
  for (nm in names(state_imm)) {
    nums <- as.integer(unlist(regmatches(nm, gregexpr("\\d+", nm))))
    ages <- if (grepl("^86", nm)) (nums[1]):max_age  # open-ended oldest group
    else if (length(nums) == 1) nums[1]       # single age (13_months, 1_year, etc.)
    else nums[1]:nums[2]                       # age range
    state_imm_exp <- c(state_imm_exp, rep(state_imm[[nm]], length(ages)))
  }
  
  state_imm_exp
  
  # expand age groups to max age
  state_imm_exp <- c()
  for (nm in names(state_imm)) {
    nums <- as.integer(unlist(regmatches(nm, gregexpr("\\d+", nm))))
    ages <- if (length(nums) == 1) (nums[1]+1):(max_age) else nums[1]:nums[2]
    state_imm_exp <- c(state_imm_exp, rep(state_imm[[nm]], length(ages)))
  }
  state_imm_exp

  # convert to ratio
  immunity_profile <- state_imm_exp/100
  
  # adjust infant immunity
  # children 6mons+ are eligible for 1 dose vaccine; 
  # children < 6mos are expected to have immunity from mother
  immunity_profile[1] <- 1/2
  
  # get suscetibility =  1 - immunity
  susceptiblilty_profile <- 1-immunity_profile
  
  # explore
  plot(susceptiblilty_profile,ylim=0:1,type='l',lwd=7,ylab='susceptibility',xlab='age')
  plot(immunity_profile,ylim=0:1,type='l',lwd=7,ylab='immunity',xlab='age')
}


# Helper to parse column names like "X0.4", "X10.14", "X85." into age ranges
parse_age_cols <- function(col_names) {
  age_cols <- col_names[grepl("^X\\d", col_names)]
  lapply(age_cols, function(nm) {
    nums <- as.integer(unlist(regmatches(nm, gregexpr("\\d+", nm))))
    if (length(nums) == 1) {
      # Open-ended like "X85." — expand to e.g. 85:89 (or adjust ceiling as needed)
      list(col = nm, ages = nums[1]:(nums[1] + 4))
    } else {
      list(col = nm, ages = nums[1]:nums[2])
    }
  })
}

expand_age_bands <- function(df) {
  non_age_cols <- df %>% select(!matches("^X\\d"))
  age_col_names <- names(df)[grepl("^X\\d", names(df))]
  
  parsed <- parse_age_cols(age_col_names)
  
  expanded <- lapply(parsed, function(p) {
    vals <- df[[p$col]]
    new_cols <- setNames(
      as.data.frame(matrix(rep(vals, length(p$ages)), ncol = length(p$ages))),
      paste0("X", p$ages)
    )
    new_cols
  })
  
  bind_cols(non_age_cols, do.call(bind_cols, expanded))
}

df_expanded <- expand_age_bands(us_imm)
