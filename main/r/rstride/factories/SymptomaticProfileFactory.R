############################################################################ #
#  This file is part of the Stride software. 
#
#  Copyright 2020, Willem L
############################################################################ #
#
# TO OBTAIN AN AGE-SPECIFIC PROPORTION OF SYMPTOMATIC CASES  
#
# APPROACH:
# 1. Use age-specific data on relative susceptibility to symptomatic infection
# 2. Impute data for children (not present)
# 3. Truncate values at 2.25 (to prevent final proportions > 1)
# 4. Set an average proportion symptomatic for the population
# 5. Rescale the age-specific profile to proportion symptomatic:
#     rescaled data =  original data / weighted mean * population average
# 6. Explore and save results
#
############################################################################ #

create_symptomatic_profile <- function()
{
        # load packages
        library(wpp2019)                                      # to obtain population data
        suppressPackageStartupMessages(library(simid.rtools)) # to save a list as XML
        
        # set overall fraction symptomatic to 50% (of something else...)
        reference_proportion <- 0.5
        
        ########################################### #
        ## LOAD AND RESHAPE DATA                 ####
        ########################################### #
        
        # Relative susceptibility to symptomatic infection by age ( Wu et al, Nature Medicine, 2020, figure 2b)
        age_profile         <- c(rep(0.16,10),  # 0-9y
                                 rep(0.16,10),  # 10-19y
                                 rep(0.4,10),   # 20-29y
                                 rep(1.0,10),   # 30-39y
                                 rep(1.3,10),   # 40-49y
                                 rep(2.0,10),   # 50-59y
                                 rep(2.9,10),   # 60-69y
                                 rep(3.0,10),   # 70-79y
                                 rep(2.5,10))   # 80y-90y
        
        # apply cutoff (to prevent final proporations > 1)
        age_profile_original          <- age_profile       # make a copy of the original data
        age_profile[age_profile>2.25] <- 2.25              # apply cutoff
        
        #' \newpage
        ########################################### #
        ## POPULATION DATA                       ####
        ########################################### #
        
        # load World Population Prospects (wpp2019) data
        data(popM)
        data_pop <- (popM[popM$name=="Belgium",c('age','2020')])
        data_pop$age_num <- seq(0,100,5)
        
        # extrapolate to all ages
        data_pop_all <- approx(x=data_pop$age_num,data_pop$`2020`,0:89,method = 'constant')
        
        # get relative proporation by age
        data_pop_rel_prop <- data_pop_all$y /sum(data_pop_all$y) * length(data_pop_all$y)
        
        ########################################### #
        ## RESCALE AGE-SPECIFIC PROPORTIONS      ####
        ########################################### #
        
        # rescale age-specific fractions
        symptomatic_profile <- age_profile / mean(age_profile*data_pop_rel_prop) * reference_proportion
        mean(symptomatic_profile)
        mean(symptomatic_profile*data_pop_rel_prop)
        
        ########################################### #
        ## EXPLORE DATA AND RESULTS              ####
        ########################################### #
        
        # open pdf
        pdf(file='./sim_output/symptomatic_profile.pdf',6,6)
        
        par(mar=c(5,5,5,1))
        
        # plot Whu data + adjustments
        plot(0:(length(age_profile_original)-1),age_profile_original,
             pch=15,
             ylab='Relative susceptibility\nto symptomatic infection',
             xlab='Age (years)',
             #main='Based on Figure 2b from\nWu et al, Nature Medicine, 2020',
             ylim=c(0,3.2)
        )
        grid(col=8)
        points(0:19,age_profile[1:20],col=2,pch=15)
        points(60:89,age_profile[61:90],col=3,pch=15)
        legend('topleft',
               c('Wu et al (2020)',
                 'Assumption (no data)',
                 'Truncated'),
               pch=15,
               col =1:3,
               bg='white',
               # title='Source',
               cex=0.8)
        
        # plot final proportions
        plot(0:(length(symptomatic_profile)-1),
             symptomatic_profile,
             ylim=c(0,1.1),type='p',lwd=3,pch=15,
             ylab='Probability symptomatic',
             xlab='Age (years)',
             #main='Probability to be symptomatic',
             xlim=c(0,90))
        #abline(h=mean(symptomatic_profile*data_pop_rel_prop),col=4,lwd=3)
        
        #abline(h=0:1)
        grid(col=8)
        # legend('topleft',
        #        'Mean (weighted)',
        #        col=4,
        #        lwd=2,
        #        bg='white',
        #        cex=0.8)
        
        # add numeric values
        text(x = which(diff(c(1,symptomatic_profile)) != 0),
             y = round(unique(symptomatic_profile),digits=2),
             round(unique(symptomatic_profile),digits=2),
             pos=3,col=1,cex=0.8)
        
        # add mean value
        # text(x = 0,
        #      y = round(mean(symptomatic_profile*data_pop_rel_prop),digits=2),
        #      round(mean(symptomatic_profile*data_pop_rel_prop),digits=2),
        #      pos=3,col=4,cex=0.8)
        
        
        dev.off()
        ########################################### #
        ## SAVE AS XML AND CSV	 	         ####
        ########################################### #
        
        # extend to 110 years of age by repeating the last one (from 0 to 110 year = 111 categories)
        max_age      <- 110
        num_age      <- max_age+1
        num_age_data <- length(symptomatic_profile)
        symptomatic_profile <- c(symptomatic_profile,rep(symptomatic_profile[num_age_data],num_age-num_age_data))
        
        # add age group as column names
        names(symptomatic_profile) <- paste0('age',0:max_age)
        
        # add info on data source and manipulation
        out_data <- unlist(list(data_source = paste('Wu et al, Nature Medicine (2020) + with assumption for children (no data) and truncated values for +60 year olds'),
                                reference_proportion = paste0('Population average: ',reference_proportion*100,'% (weighted by age)'),
                                data_processing    = "Rescaled age-specific relative susceptibility to symptomatic infection to proportion symptomatic",
                                round(symptomatic_profile,digits=4)))
        
        # save as xml (with meta data)
        smd_save_as_xml(data_list = out_data, root_name = 'prob_symptomatic', file_name_prefix = 'symptomatic_covid19',
                        xml_prefix = paste0(' This file is part of the Stride software [', format(Sys.time()), ']'))
        
        # save as csv (only the values)
        write.table(data.frame(age = 0:max_age,proportion = symptomatic_profile,population_average = reference_proportion),
                    file = 'symptomatic_covid19.csv',sep=',',row.names=F)
        
} # end function

