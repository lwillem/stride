#############################################################################
#
# PREPARE USA-SPECIFIC POPULATION FILE FOR STRIDE
# 
############################################################################ #

# to debug/develop
if(0==1){
  dane <- getFREDdata(state = "WI", county = "Dane")
  mke  <- getFREDdata(state = "WI", county = "Milwaukee")
  cook <- getFREDdata(state = "IL", county = "Cook")
  hennepin <- getFREDdata(state = "MN", county = "Hennepin")
  
  # or set function parameters
  state = "WI"; county = "Dane"
}

options(scipen=999)    # turn off scientific notation
# options(scipen=0)      # turn on scientific notation

############################### Synthetic FRED Population ################################

## Manually download data from FRED: https://fred.publichealth.pitt.edu/syn_pops to '~/opt/FRED_population_usa/'

# MN - Hennepin
# MKE - Milwaukee
# IL - Cook

getFREDdata <- function(state = "WI", county = "Milwaukee", com_target_size = 1000, rng_seed = 1234){
  
  # set rng seed
  set.seed(rng_seed)
  
  print("FRED population files for the desired county must be manually downloaded from 'https://fred.publichealth.pitt.edu/syn_pops'
  and stored in the '~/opt/FRED_population_usa/' folder before running this script.")
  
  folder <- fips(state, county = county)
  state_name_full <- fips_info(folder)$full
  cty_name_full <- fips_info(folder)$county
  
  dir <- grep(folder, list.dirs(paste0("~/opt/FRED_population_usa/", state)), value = TRUE)

  if(length(dir) == 0){
    stop("!!! County data not found !!! \nCheck if population files are stored in '~/opt/FRED_population_usa' 
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
                             com_target_size      = com_target_size,     # target size of the community
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
  
  return(pop_data)
}

# Calculate the distance for two [X,Y] coordinates
get_distance <- function(f_x_coord1,f_y_coord1,f_x_coord2,f_y_coord2){
  sqrt((f_x_coord1-f_x_coord2)^2 + (f_y_coord1-f_y_coord2)^2)
  
}

# Assign community pools
# Derived from the STRIDE POPULATION repository
generate_population_community <- function(pop_data, pop_settings, district_data, district_distance)
{
  
  # INITIALISE
  pop_data$com_id  <- 0
  pop_data$com2_id <- 0
  
  ## SET TARGET SIZE
  com_target_size       <- pop_settings$com_target_size
  num_com               <- pop_settings$size / com_target_size
  if(com_target_size == 0){ num_com <- 0}
  
  ## SET COMMUNITY LOCATIONS BASED ON POPULATION DENSITY
  # sample people as "center of a community"
  sample_population     <- sample(pop_settings$size,num_com,replace=F)
  community_data        <- pop_data[sample_population,c('stcotrbg','latitude','longitude')]
  
  if(num_com > 0){  #special case
    
    pop_data$com_distance <- NA
    pop_data$com2_distance <- NA
    # ASSIGN PEOPLE TO COMMUNITIES
    
    # PRE-GROUP: build district -> row index list ONCE (avoids repeated full-vector scans)
    district_idx <- split(seq_len(nrow(pop_data)), pop_data$stcotrbg)
    
    # PRE-INDEX: lookup district coords by GEOID10 ONCE
    district_lookup <- setNames(seq_len(nrow(district_data)), district_data$GEOID10)
    
    for(i_district in seq_len(nrow(district_distance))){
      bg10 <- rownames(district_distance)[i_district]
      
      idx <- district_idx[[bg10]]
      if(is.null(idx)) next
      district_size <- length(idx)
      if(district_size == 0) next
      
      dd_row <- district_lookup[bg10]
      district_x_coord <- district_data$latitude[dd_row]
      district_y_coord <- district_data$longitude[dd_row]
      
      distances <- get_distance(district_x_coord, district_y_coord,
                                community_data$latitude, community_data$longitude)
      
      sel_com_id   <- order(distances)[1:3]
      sel_com_prob <- rep(1, length(sel_com_id))
      sel_com_prob <- sel_com_prob / sum(sel_com_prob)
      
      pop_data$com_id[idx]  <- sample(sel_com_id, district_size, replace = TRUE, prob = sel_com_prob)
      pop_data$com2_id[idx] <- sample(sel_com_id, district_size, replace = TRUE, prob = sel_com_prob)
      
      pop_data$com_distance[idx]  <- distances[pop_data$com_id[idx]]
      pop_data$com2_distance[idx] <- distances[pop_data$com2_id[idx]]
    }
    
    # names(pop_data)
    ## UPDATE COMMUNITY DATA
    community_data$size1 <- 0
    community_data$size2 <- 0 
    tmp <- table(pop_data$com_id)
    community_data$size1[as.numeric(names(tmp))] <- tmp 
    tmp <- table(pop_data$com2_id)
    community_data$size2[as.numeric(names(tmp))] <- tmp 
    
    plot(sort(table(pop_data$com_id)),xlab='id',ylab='size',lwd=6)
    lines(sort(table(pop_data$com2_id)),col=2)
    abline(h=com_target_size,col=4)
    legend('topright',c('community 1', 'community 2', 'target size'),col=c(1,2,4),lwd=2)
    
    # adapt community_data
    community_data <- data.frame(id=1:nrow(community_data),
                                 district_id=community_data$stcotrbg, 
                                 size1=community_data$size1, 
                                 size2=community_data$size2)
    names(community_data)
    
    # add x and y coordinates
    community_data <- merge(community_data, 
                            data.frame(district_id = district_data$GEOID10,
                                       x_coord = district_data$latitude,
                                       y_coord = district_data$longitude))
    
    
    ## EXPLORE DISTANCE FROM HOME
    tbl_dist <- table(cut(c(pop_data$com_distance, pop_data$com2_distance), c(0,.002,.004,.006,.02), right = FALSE))
    barplot(tbl_dist / sum(tbl_dist), 
            xlab = 'Distance from home', 
            ylab = 'Fraction of the population')
    
    # FAST community overlap: single crosstab instead of per-community loop
    overlap_tbl <- table(factor(pop_data$com_id, levels = 1:num_com),
                         factor(pop_data$com2_id, levels = 1:num_com))
    opt <- matrix(overlap_tbl, nrow = num_com, ncol = num_com)
    
    plot(sort(opt), main = paste0('community overlap (', num_com, 'x', num_com, ')'),
         xlab = "Community (sorted)", ylab = 'Overlap: count')
    hist(opt[opt > 0], 500, xlab = 'Community overlap (if >0)', main = 'Community overlap (if >0)')
    
  } # end if community size > 0
  
  ## EXPLORE
  df_sf <- st_as_sf(district_data, coords = c("longitude", "latitude"), crs = 4326)
  community_data_sf <- st_as_sf(community_data, coords = c("y_coord", "x_coord"), crs = 4326)
  gmap <- ggplot() +
    geom_sf(data = df_sf, color = "red", size = 2) +
    geom_sf(data = community_data_sf, color = "blue", size = 1) +
    theme_minimal()
  print(gmap)
  
  barplot(table(cut(community_data$size1, seq(0,2,0.2) * com_target_size)) / num_com,
          xlab = 'community size', ylab = 'fraction', las = 2, cex.names = 0.6)
  
  barplot(rbind(sort(table(pop_data$com_id)), sort(table(pop_data$com2_id))),
          beside = TRUE, xlab = 'community (sorted)', ylab = 'size', xaxt = 'n')
  abline(h = com_target_size)
  legend('topleft', c('community 1','community 2'), fill = c(1,8))
  
  hist(c(table(pop_data$com_id), table(pop_data$com2_id)), xlab = 'community size')
  
  # return pop_data and community_data
  return(list(pop_data,community_data))
  
}

