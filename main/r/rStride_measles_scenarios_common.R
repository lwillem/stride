#############################################################################
#  This file is part of the Stride software. 
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
#
#  Copyright 2026, Manansala R, Willem L
#############################################################################
#
# Helpers for the county measles scenario scripts (rStride_measles_<scenario>.R).
# No scenario or policy parameters are set here: those are in each scenario script.
# This file only holds the county input files and the loop that runs one rSTRIDE
# project per county.
#
#############################################################################

# Load rStride
source('./bin/rstride/rStride.R')

# county input files (the same in every scenario)
get_measles_county_files <- function(county){
  files <- list(
    CT = c("pop_usa_ct_fairfield_c1000_hhdense.csv",   "contact_matrix_usa_ct_fairfield_c1000.xml",   "CT_20261022_164142_r0_disease_measles_usa.xml"),
    ME = c("pop_usa_me_piscataquis_c1000_hhdense.csv", "contact_matrix_usa_me_piscataquis_c1000.xml", "ME_20261023_081954_r0_disease_measles_usa.xml"),
    MN = c("pop_usa_mn_hennepin_c1000_hhdense.csv",    "contact_matrix_usa_mn_hennepin_c1000.xml",    "MN_20261022_192758_r0_disease_measles_usa.xml"),
    OH = c("pop_usa_oh_cuyahoga_c1000_hhdense.csv",    "contact_matrix_usa_oh_cuyahoga_c1000.xml",    "OH_20261022_232515_r0_disease_measles_usa.xml"),
    TX = c("pop_usa_tx_gaines_c1000_hhdense.csv",      "contact_matrix_usa_tx_gaines_c1000.xml",      "TX_20261023_111153_r0_disease_measles_usa.xml"),
    SC = c("pop_usa_sc_spartanburg_c1000_hhdense.csv", "contact_matrix_usa_sc_spartanburg_c1000.xml", "SC_20261019_191754_r0_disease_measles_usa.xml")
  )
  f <- files[[county]]
  if(is.null(f)) stop("unknown county: ", county)
  list(population_file            = file.path("data/population_6region", f[1]),
       age_contact_matrix_file    = file.path("data/contacts_usa", f[2]),
       disease_config_file        = file.path("data", f[3]),
       immunity_distribution_file = sprintf("data/immunity_6region/immunity_measles_adult_%s_calib_f035.xml", county),
       vaccine_distribution_file  = sprintf("data/immunity_6region/immunity_measles_child_%s.xml", county))
}

# Run one rSTRIDE project per row of county_param.
# exp_param_list: parameters shared by all counties
# county_param:   data.frame with a column 'county' and one column per parameter to set for that county
#                 (column names = rSTRIDE parameter names)
run_measles_counties <- function(exp_param_list, county_param, scenario_name,
                                 num_parallel_workers = 3, inspect_output = TRUE){

  project_dirs <- list()
  for(i in seq_len(nrow(county_param))){
    county <- county_param$county[i]

    county_param_list <- exp_param_list
    county_files      <- get_measles_county_files(county)
    for(p in names(county_files))                 county_param_list[[p]] <- county_files[[p]]
    for(p in setdiff(names(county_param),'county')) county_param_list[[p]] <- county_param[[p]][i]
    if(county_param_list$vaccine_hesitancy_num_clusters == 0){
      county_param_list$vaccine_hesitancy_cluster_pooltype <- "None"
    }

    # get grid-based design of experiments
    exp_design <- .rstride$get_full_grid_exp_design(exp_param_list = county_param_list,
                                                    num_rng_seeds  = county_param_list$num_rng_seeds)

    # run rSTRIDE (remove_run_output must stay FALSE to get the household aggregation)
    project_dir <- run_rStride(exp_design,
                               paste0('_measles_', county, '_', scenario_name),
                               remove_run_output    = FALSE,
                               num_parallel_workers = num_parallel_workers)

    if(inspect_output){
      inspect_summary(project_dir)
      inspect_participant_data(project_dir)
      inspect_transmission_dynamics(project_dir)
      inspect_incidence_data(project_dir)
      inspect_prevalence_data(project_dir)
    }
    project_dirs[[county]] <- project_dir
  }
  return(project_dirs)
}
