#############################################################################
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
#  Copyright 2018, Willem L
#############################################################################
#
# PREPARE COMMUNITY DATA - ADAPTED FROM STRIDE
# For simulated household data from FRED
# using US PUMS and TIGRIS Data (R package 'tigris')
#
#############################################################################

# Calculate the distance for two [X,Y] coordinates
get_distance <- function(f_x_coord1,f_y_coord1,f_x_coord2,f_y_coord2){
  sqrt((f_x_coord1-f_x_coord2)^2 + (f_y_coord1-f_y_coord2)^2)
  
}

generate_population_community <- function(pop_data,pop_settings,district_data)
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
    # per district, select one of the 3 nearest communities
    i_district <- 2
    for(i_district in 1:nrow(district_distance)){
      # get people from one district and temporary store district features
      bg10 <- rownames(district_distance)[i_district]
      flag_district    <- pop_data$stcotrbg == bg10
      district_x_coord <- district_data[district_data$GEOID10 == bg10, "latitude"] #unique(pop_data$latitude[flag_district])
      district_y_coord <- district_data[district_data$GEOID10 == bg10, "longitude"] #unique(pop_data$longitude[flag_district])
      district_size    <- sum(flag_district)
      
      if(district_size > 0){
        # get distances from district 'i' to all communities
        distances <- get_distance(district_x_coord$latitude,district_y_coord$longitude,community_data$latitude,community_data$longitude)
        # select com_id from nearest 3 communities
        sel_com_id <- which(distances %in% (sort(unique(distances))[1:3]))
        # calculate probability based on distance ranking: 1/rank
        sel_com_prob <- rep(1,length(sel_com_id))
        # sel_com_prob <- 1/order(distances[sel_com_id])
        # standardize probabilities
        sel_com_prob <- sel_com_prob / sum(sel_com_prob)
        # sample from community selection, with replacement
        pop_data$com_id[flag_district]  <- sample(sel_com_id,district_size,replace=T,prob=sel_com_prob)
        pop_data$com2_id[flag_district] <- sample(sel_com_id,district_size,replace=T,prob=sel_com_prob)
        
        pop_data$com_distance[flag_district] <- distances[pop_data$com_id[flag_district]]
        pop_data$com2_distance[flag_district] <- distances[pop_data$com2_id[flag_district]]
        
      }
    }
    
    names(pop_data)
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
    tbl_dist <- table(cut(c(pop_data$com_distance,pop_data$com2_distance),c(0,.002,.004,.006,.02),right=F))
    barplot(tbl_dist/sum(tbl_dist),
            xlab='Distance from home',
            ylab='Fraction of the population')
    
    ## EXPLORE COMMUNITY OVERLAP
    opt <- matrix(0,nrow=num_com,ncol=num_com)
    i <- 1
    for(i in 1:num_com){
      tmp <- table(pop_data$com2_id[pop_data$com_id == i])
      opt[i,as.numeric(names(tmp))] <- tmp
    }
    
    plot(sort(opt),
         main=paste0('community overlap (',num_com,'x',num_com,')'),
         xlab="Community (sorted)",
         ylab='Overlap: count')
    
    hist(opt[opt>0],500,xlab='Community overlap (if >0)',
         main='Community overlap (if >0)')
    
  } # end if community size > 0
  
  ## EXPLORE
  df_sf <- st_as_sf(district_data, coords = c("longitude", "latitude"), crs = 4326)
  community_data_sf <- st_as_sf(community_data, coords = c("y_coord", "x_coord"), crs = 4326)
  gmap <- ggplot() +
    # geom_sf(data = world, fill = "gray95") +
    geom_sf(data = df_sf, color = "red", size = 2) +
    geom_sf(data = community_data_sf, color = "blue", size = 1) +
    theme_minimal()
  print(gmap)
  
  barplot(table(cut(community_data$size1,seq(0,2,0.2)*com_target_size))/num_com,
          xlab='community size',
          ylab='fraction',
          las=2,
          cex.names = 0.6)
  
  barplot(rbind(sort(table(pop_data$com_id)),
                sort(table(pop_data$com2_id))),
          beside=T,
          xlab='community (sorted)',
          ylab='size',
          xaxt='n')
  abline(h=com_target_size)
  legend('topleft',c('community 1','community 2'),fill=c(1,8))
  
  hist(c((table(pop_data$com_id)),
         (table(pop_data$com2_id))),
       xlab='community size')
  
  # return pop_data and community_data
  return(list(pop_data,community_data))
  
}
