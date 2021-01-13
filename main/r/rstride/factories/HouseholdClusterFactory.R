############################################################################ #
#  This file is part of the Stride software. 
#
#  Copyright 2020, Willem L
############################################################################ #
#
# TO OBTAIN AGE-SPECIFIC HOUSEHOLD CLUSTERS  
#
# APPROACH:
# 1. Use pre-computed populaitn data
# 2. Match households by the age the oldest member
# 3. Save new population data
#
############################################################################ #


if(0==1){

    # set filename
    pop_file_name <- 'pop_belgium600k_c500_teachers_censushh.csv'
    
    # set maximum age difference between household seniors.
    max_age_diff <- 3
    
    # max number of households in one cluster?
    household_cluster_size <- 2
    
}

extend_population_data <- function(pop_file_name,max_age_diff,household_cluster_size){

smd_print("START HOUSEHOLD CLUSTERING")
        
########################################### #
## LOAD DATA ----
########################################### #

# load data
pop_data <- read.table(smd_file_path('./data',pop_file_name),sep=',',header=T)

# inspection
names(pop_data)
names_orig <- names(pop_data)

# get file names for output files
pop_file_name_out <- paste0('sim_output/',gsub('.csv',paste0('_extended',max_age_diff,'_size',household_cluster_size,'.csv'),pop_file_name))
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
    smd_print("CLUSTER SIZE TO SMALL... STOP")
    return(NULL)
}

if(num_other_households == 0){
    # create new variable for household_cluster_id
    pop_data$household_cluster_id <- pop_data$household_id
       
} else {
    
    # get info on max age per hh
    hh_data_age <- aggregate(age ~ household_id ,data=pop_data,max)     # aggregate
    dim(hh_data_age)
    
    # add primary community
    # note: households are duplicated if two seniors have the same age (and different prim. community)
    hh_data_summary <- merge(hh_data_age,pop_data[,c('primary_community',names(hh_data_age))])
    dim(hh_data_summary)
    
    # select one row per household
    # sort by household id, and remove duplicated rows
    hh_data_summary <- hh_data_summary[order(hh_data_summary$household_id),]
    flag_unique     <- hh_data_summary$household_id != c(hh_data_summary$household_id[-1],0)
    table(flag_unique)
    hh_data_summary <- hh_data_summary[flag_unique,]
    
    # select unique primary community ids
    community_opt <- sort(unique(hh_data_summary$primary_community))
    length(community_opt)
    
    # create new variable for household_cluster_id
    hh_data_summary$household_cluster_id <- NA
    
    # start parallel nodes
    par_nodes_info <- smd_start_cluster()
    
    tmp_time <- Sys.time()
    i_community <- community_opt[1]
    # add household cluster, based on the primary community of each household senior
    foreach(i_community = community_opt,
            .export = c('smd_print_progress',rStride_functions),
            .combine = rbind) %dopar%
    {
            # print progress        
            smd_print_progress(i_community,length(community_opt),tmp_time,par_nodes_info)        
    
            # select population data
            hh_data_community <- hh_data_summary[hh_data_summary$primary_community == i_community,]
            dim(hh_data_community)
            
            # start with cluster id 1
            hh_cluster_counter <- 1
            
            i_hh <- 1
            # loop over the households... and match
            for(i_hh in 1:nrow(hh_data_community)){
                    print(paste(i_hh,hh_cluster_counter))
                    
                    if(is.na(hh_data_community$household_cluster_id[i_hh])){
                            
                            household_id   <- hh_data_community$household_id[i_hh]
                            hh_age_max     <- hh_data_community$age[i_hh]
                            
                            flag_opt <- hh_data_community$age %in% seq(hh_age_max-max_age_diff,hh_age_max+max_age_diff) &
                                        hh_data_community$household_id != household_id &
                                        is.na(hh_data_community$household_cluster_id)
                            table(flag_opt)
                            
                            if(any(flag_opt)){
                                
                                # set cluster size and sample
                                c_size <- min(sum(flag_opt),household_cluster_size-1)
                                household_id_sample <- sample(hh_data_community$household_id[flag_opt],c_size)
                                
                                # select households
                                flag_cluster <- hh_data_community$household_id %in% c(household_id,household_id_sample)
                                hh_data_community$household_cluster_id[flag_cluster] <- paste0(i_community,'-',hh_cluster_counter)
                                        
                                # update
                                hh_cluster_counter = hh_cluster_counter + 1
                            } # end if-clause: any in flag_opt
                    } # end if-clause: household is not part a cluster (yet)
            } # end for-loop: households within one community
            # return
            return(hh_data_community)
    } -> hh_data_summary # end parallel foreach
    
    # close parallel nodes
    smd_stop_cluster()
    
    # reformat household_cluster_id into numeric value
    hh_data_summary$household_cluster_id <- as.numeric(as.factor(hh_data_summary$household_cluster_id))
    hh_data_summary$household_cluster_id[is.na(hh_data_summary$household_cluster_id )] <- 0
    
    # add household_cluster_id to pop_data
    names_pop_data   <- names(pop_data)
    names_hh_summary <- c('household_id','household_cluster_id')
    pop_data <- merge(pop_data,hh_data_summary[,names_hh_summary])[, union(names_pop_data, names_hh_summary)]
    names(pop_data)
    
    # check pop_data
    head(pop_data)
    table(is.na(pop_data$household_cluster_id))
    table(pop_data$household_cluster_id == 0) / nrow(pop_data)
    
    # check cluster size
    table(table(pop_data$household_cluster_id))
    hist(table(pop_data$household_cluster_id[pop_data$household_cluster_id!= 0]))
    
    
    # example
    pop_data[pop_data$household_cluster_id == '75689',]
} # end if-else cluster-size is 1

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


