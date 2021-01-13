############################################################################ #
# PARALLEL ADMIN
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Copyright (C) 2020 lwillem, SIMID, UNIVERSITY OF ANTWERP, BELGIUM
############################################################################ #
#' @import doParallel
#' @import parallel
#' @import foreach
#' @import iterators


############################################################### #
## START CLUSTER WITH PARALLEL WORKERS ----
############################################################### #
#' @title Start parallel working nodes
#' @param timeout The timeout in seconds for the nodes
#' @param num_proc The number of parallel workers to start (NA = default, number of CPU cores)
#' @return Cluster info, also available as global variable 'par_nodes_info'
#' @keywords external
#' @export
#' @note Explicitly setup the PSOCK cluster using the 'sequential' setup strategy, to prevent an error in RStudio on macOS with R 4.0.x (https://github.com/rstudio/rstudio/issues/6692)
smd_start_cluster <- function(timeout = 100, num_proc = NA)
{
  smd_print("START PARALLEL WORKERS")

  if(is.na(num_proc)){
    num_proc <- detectCores()
  }

  # check num_proc max and min... and use default if needed
  if(num_proc > detectCores() || num_proc < 1){
    num_proc <- detectCores()
    smd_print("Number of workers given to 'smd_start_cluster()' is not valid, so used number of CPU cores",WARNING = T)
  }

  ## SETUP PARALLEL NODES
  # note: they will be removed after 100 seconds inactivity
  par_cluster   <- makeCluster(num_proc, cores=num_proc, timeout = timeout, setup_strategy = "sequential")
  registerDoParallel(par_cluster)

  # store the process id (pid) of the first slave
  pid_slave1 <- clusterEvalQ(par_cluster, { Sys.getpid() })[[1]]

  # CREATE GLOBAL VARIABLE
  par_nodes_info <<- list(par_cluster = par_cluster,
                          pid_master  = Sys.getpid(),
                          pid_slave1  = pid_slave1)
}


############################################################### #
## STOP CLUSTER WITH PARALLEL WORKERS ----
############################################################### #
#' @title Stop parallel working nodes
#' @keywords external
#' @export
smd_stop_cluster <- function()
{
  ## CLOSE NODES AND NODE INFO
  if(exists('par_nodes_info')){
    smd_print("STOP PARALLEL WORKERS")

    stopCluster(par_nodes_info$par_cluster);
    rm(par_nodes_info, envir = .GlobalEnv) # REMOVE GLOBAL VARIABLE
  }
}

############################################################### #
## CHECK IF CLUSTER EXISTS AND START ONE IF NOT PRESENT----
############################################################### #
#' @title Check parallel working nodes, and start a cluster if not present yet.
#' @return Cluster info
#' @export
smd_check_cluster <- function()
{
  if(!exists('par_nodes_info')){
    smd_start_cluster()
  } else if (!any(grepl(par_nodes_info$pid_slave1,system('ps -A',intern = T)))){
    smd_start_cluster()
  }

  # return par_nodes_info
  return(par_nodes_info)
}
