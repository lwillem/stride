############################################################################ #
#  This file is part of the Stride software. 
#
#  Copyright 2024, Willem L
############################################################################ #
#
# TO OBTAIN AGE-SPECIFIC HOUSEHOLD CLUSTERS  
#
# APPROACH:
# 1. Use pre-computed population data
# 2. Match households by the age the oldest member
# 3. Save new population data
#
############################################################################ #

# load simid.rtools package 
#suppressPackageStartupMessages(library(simid.rtools)) # to save a list as XML

if(0==1){
  
  # set filename
  pop_file_name <- '20260825_094649_population_TX-Gaines.csv'
  
  # set maximum age difference between household seniors.
  max_age_diff <- 3
  
  # max number of households in one cluster?
  household_cluster_size <- 10
  
  # percent clustering
  # pct_households_clustered <- 0.0499   # CT
  # pct_households_clustered <- 0.105    # ME
  # pct_households_clustered <- 0.0619     # MN
  # pct_households_clustered <- 0.1219     # OH
  # pct_households_clustered <- 0.1443     # SC
  pct_households_clustered <- 0.1166     # TX

  # set seed
  seed <- 1234567
  
  # call function
  extend_population_data(pop_file_name            = pop_file_name,
                         max_age_diff             = max_age_diff,
                         household_cluster_size   = household_cluster_size,
                         pct_households_clustered = pct_households_clustered,
                         seed                     = seed)
  
}

extend_population_data <- function(pop_file_name, max_age_diff, household_cluster_size,
                                   pct_households_clustered, seed = NA){
  
  smd_print("START HOUSEHOLD CLUSTERING")
  
  ########################################### #
  ## LOAD DATA ----
  ########################################### #
  
  # load data
  pop_data <- read.table(smd_file_path('./data/population_6region/',pop_file_name),
                         sep=',',header=T)
  
  # inspection
  names(pop_data)
  names_orig <- names(pop_data)
  
  # get file names for output files
  pop_file_name_out <- paste0('data/population_6region/',
                              gsub('.csv',paste0('_extended', max_age_diff,
                                                 '_size', household_cluster_size,
                                                 '_pct', round(pct_households_clustered*100,0),
                                                 '.csv'),pop_file_name))
  pop_file_zip_out  <- gsub('.csv','.zip',pop_file_name_out)
  
  # CLI statement
  smd_print(pop_file_name_out)
  
  ########################################### #
  ##  HOUSEHOLD EXTENSION ----
  ########################################### #
  # if every HH teams up with other HH(s)?
  
  # define the total cluster size 
  num_other_households <- household_cluster_size-1
  
  if(num_other_households < 0){
    smd_print("CLUSTER SIZE TOO SMALL... STOP")
    return(NULL)
  }
  
  if(num_other_households == 0){
    # create new variable for household_cluster_id
    pop_data$household_cluster_id <- pop_data$household_id
    
  } else {
    
    # get info on max age per hh
    hh_data_age <- aggregate(age ~ household_id, data=pop_data, max)   # aggregate
    dim(hh_data_age)
    
    # add primary community
    # note: households are duplicated if two seniors have the same age (and different prim. community)
    hh_data_summary <- merge(hh_data_age, pop_data[,c('primary_community', names(hh_data_age))])
    dim(hh_data_summary)
    
    # select one row per household
    # sort by household id, and remove duplicated rows
    hh_data_summary <- hh_data_summary[order(hh_data_summary$household_id),]
    flag_unique     <- hh_data_summary$household_id != c(hh_data_summary$household_id[-1],0)
    table(flag_unique)
    hh_data_summary <- hh_data_summary[flag_unique,]
    
    num_households <- nrow(hh_data_summary)
    
    # NEW: compute the target household count from the percentage and this
    # population's actual household count, so pct_households_clustered can
    # be reused as-is across population files of differing size
    target_households_clustered <- round(num_households * pct_households_clustered)
    
    if(target_households_clustered > num_households){
      smd_print("TARGET EXCEEDS TOTAL HOUSEHOLDS... STOP")
      return(NULL)
    }
    
    if(!is.na(seed)) set.seed(seed)
    
    smd_print(paste("Target households to cluster:", target_households_clustered,
                    "out of", num_households,
                    sprintf("(%.1f%%)", pct_households_clustered*100)))
    
    # select unique primary community ids, shuffled order
    community_opt <- sample(unique(hh_data_summary$primary_community))
    length(community_opt)
    
    # create new variable for household_cluster_id
    hh_data_summary$household_cluster_id <- NA
    
    tmp_time <- Sys.time()
    total_clustered <- 0
    cluster_counter <- 0
    num_skipped_no_match <- 0   # households that couldn't find an age-window match
    
    # add household cluster, based on the primary community of each household senior
    for(i_community in community_opt){
      
      if(total_clustered >= target_households_clustered) break
      
      # select population data
      hh_data_community <- hh_data_summary[hh_data_summary$primary_community == i_community,]
      dim(hh_data_community)
      
      # sort this community's households by age (oldest hh member)
      ord         <- order(hh_data_community$age)
      idx_sorted  <- hh_data_community[ord, "household_id"]
      ages_sorted <- hh_data_community[ord, "age"]
      n_hh        <- length(idx_sorted)
      
      i_hh <- 1 # counter to iterate over num hh in a community (n_hh)
      # loop over the households... and match
      while(i_hh <= n_hh && total_clustered < target_households_clustered){
        
        remaining_needed <- target_households_clustered - total_clustered
        
        # get number of households to cluster in community 
        # if remaining hh not in community < household_cluster_size, then find leftover households 
        # or get number of household clusters needed to hit target
        max_chunk_size <- min(household_cluster_size, n_hh - i_hh + 1, remaining_needed)
        
        # NEW: shrink chunk to respect max_age_diff. Since the list is sorted
        # by age, the age range of a chunk is simply (last age - first age)
        chunk_size <- max_chunk_size
        while(chunk_size >= 2 &&
              (ages_sorted[i_hh + chunk_size - 1] - ages_sorted[i_hh]) > max_age_diff){
          chunk_size <- chunk_size - 1
        }
        
        # skip forming a cluster of size 1 (no partner) unless it's genuinely
        # the last household needed to hit the target exactly
        if(chunk_size < 2){
          if(remaining_needed < 2){
            break # target basically satisfied, nothing more to add
          }
          # this household has no valid age-window partner starting here ->
          # move on to the next household in this community
          num_skipped_no_match <- num_skipped_no_match + 1
          i_hh <- i_hh + 1
          next
        }
        
        idx_chunk <- idx_sorted[i_hh:(i_hh + chunk_size - 1)]
        cluster_counter <- cluster_counter + 1
        hh_data_summary[hh_data_summary$household_id %in% idx_chunk, "household_cluster_id"] <- paste0(i_community, '-', cluster_counter)
        
        total_clustered <- total_clustered + chunk_size
        i_hh <- i_hh + chunk_size
      }
    }
    
    # sanity checks
    smd_print(sprintf("Target: %d households (%.2f%%) | Actually clustered: %d (%.2f%%)",
                      target_households_clustered, 100*pct_households_clustered,
                      total_clustered, 100*total_clustered/num_households))
    smd_print(sprintf("Households skipped due to max_age_diff constraint: %d", num_skipped_no_match))
    
    if(total_clustered < target_households_clustered){
      smd_print(sprintf("WARNING: could not reach target - ran out of eligible age-matched households. Short by %d.",
                        target_households_clustered - total_clustered))
      smd_print("Consider increasing max_age_diff, increasing household_cluster_size, or lowering pct_households_clustered.")
    }
    
    # reformat household_cluster_id into numeric value
    hh_data_summary$household_cluster_id <- as.numeric(as.factor(hh_data_summary$household_cluster_id))
    hh_data_summary$household_cluster_id[is.na(hh_data_summary$household_cluster_id )] <- 0
    
    # add household_cluster_id to pop_data
    names_pop_data   <- names(pop_data)
    names_hh_summary <- c('household_id','household_cluster_id')
    pop_data <- merge(pop_data,hh_data_summary[,names_hh_summary])[, union(names_pop_data, names_hh_summary)]
    names(pop_data)
    
    # check actual clustered % achieved
    pct_actual <- mean(pop_data$household_cluster_id != 0)
    smd_print(sprintf("Actual %% of individuals in a household cluster: %.1f%%", pct_actual*100))
    
    # check pop_data
    head(pop_data)
    table(is.na(pop_data$household_cluster_id))
    table(pop_data$household_cluster_id == 0) / nrow(pop_data)
    
    # check cluster size
    table(table(pop_data$household_cluster_id))
    hist(table(pop_data$household_cluster_id[pop_data$household_cluster_id!= 0]))
    
    # check cluster size (in households)
    hh_per_cluster <- table(hh_data_summary$household_cluster_id[hh_data_summary$household_cluster_id != 0])
    smd_print("Distribution of cluster sizes (in households):")
    print(table(hh_per_cluster))
    
    
    # example
    # pop_data[pop_data$household_cluster_id == '75689',]
  } # end if-else cluster-size is 1
  
  # Diagnostics
  cluster_diagnostics(pop_data, pop_file_name_out)
  
  # STORE CSV FILE  ####
  #################### #
  write.table(pop_data,file = pop_file_name_out,sep=',',row.names=F)
  
  # STORE ZIP ARCHIVE ####
  ###################### #
  zip(zipfile = paste0(pop_file_zip_out),
      files = pop_file_name_out,
      flags = "-r9Xj")
  
  
  smd_print("HOUSEHOLD CLUSTERING COMPLETE")
  
} # end function

##  CLUSTER AGE DIAGNOSTICS
cluster_diagnostics <- function(pop_data, pop_file_name_out){

  clustered_pop <- pop_data[pop_data$household_cluster_id != 0, ]

  if(nrow(clustered_pop) > 0){

    # per-cluster: oldest / youngest individual, num households, num members,
    # presence of children (age < 18), num households with at least one child

    cluster_age_summary <- aggregate(age ~ household_cluster_id, data = clustered_pop,
                                     FUN = function(x) c(min = min(x), max = max(x)))
    cluster_age_summary <- do.call(data.frame, cluster_age_summary)
    names(cluster_age_summary) <- c('household_cluster_id', 'youngest_member_age', 'oldest_member_age')

    # number of unique households per cluster
    hh_count <- aggregate(household_id ~ household_cluster_id, data = clustered_pop,
                          FUN = function(x) length(unique(x)))
    names(hh_count) <- c('household_cluster_id', 'num_households')

    # total individuals per cluster
    member_count <- aggregate(household_id ~ household_cluster_id, data = clustered_pop,
                              FUN = length)
    names(member_count) <- c('household_cluster_id', 'num_individuals')

    # households with at least one child (age < 18) per cluster
    clustered_pop$is_child <- clustered_pop$age < 18
    hh_child_flag <- aggregate(is_child ~ household_id + household_cluster_id, data = clustered_pop,
                               FUN = any)
    hh_with_children <- aggregate(is_child ~ household_cluster_id, data = hh_child_flag,
                                  FUN = sum)
    names(hh_with_children) <- c('household_cluster_id', 'num_households_with_children')

    # age range spread of household-senior ages within the cluster
    # (reflects how well max_age_diff constraint held; senior age = per household max age)
    hh_senior_ages <- aggregate(age ~ household_id + household_cluster_id, data = clustered_pop, FUN = max)
    senior_age_range <- aggregate(age ~ household_cluster_id, data = hh_senior_ages,
                                  FUN = function(x) max(x) - min(x))
    names(senior_age_range) <- c('household_cluster_id', 'senior_age_range')

    # combine all diagnostics into one table
    cluster_diagnostics <- Reduce(function(x,y) merge(x,y,by='household_cluster_id'),
                                  list(cluster_age_summary, hh_count, member_count,
                                       hh_with_children, senior_age_range))
    cluster_diagnostics$pct_households_with_children <- round(100 * cluster_diagnostics$num_households_with_children /
                                                                cluster_diagnostics$num_households, 1)

    smd_print("CLUSTER AGE DIAGNOSTICS (first rows):")
    print(head(cluster_diagnostics))

    smd_print("Summary across all clusters:")
    print(summary(cluster_diagnostics[, c('youngest_member_age','oldest_member_age',
                                          'num_households','num_individuals',
                                          'senior_age_range','pct_households_with_children')]))

    smd_print(sprintf("Clusters exceeding max_age_diff (%d) in senior age range: %d out of %d",
                      max_age_diff,
                      sum(cluster_diagnostics$senior_age_range > max_age_diff),
                      nrow(cluster_diagnostics)))

    # save diagnostics alongside population output
    diagnostics_file_out <- gsub('.csv', '_cluster_diagnostics.csv', pop_file_name_out)
    write.table(cluster_diagnostics, file = diagnostics_file_out, sep=',', row.names=F)
    smd_print(paste("Cluster diagnostics saved to:", diagnostics_file_out))
  }

} # end if-else cluster-size is 1
# 
