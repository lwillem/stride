#############################################################################
# FILE PATH ADMIN
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
# Copyright (C) 2019 lwillem, SIMID, UNIVERSITY OF ANTWERP, BELGIUM
#############################################################################

#############################################################################
# LOAD RDATA FILE WITH A SINGLE VARIABLE
#############################################################################

#' @title Reads a file in RData format and returns the stored variable.
#' @param filename    The filename of the .RData File (case sensitive!)
#' @keywords external
#' @export
#' @note This function enables to rename a stored variable in the (compressed)
#' RData format. If the file contains multiple variables, this functions returns NULL.
smd_read_rdata <- function(filename){

    # check file extension
  filename_extension <- tolower(substr(filename,nchar(filename)-5,nchar(filename)))
  if(filename_extension != '.rdata'){
    smd_print('The given file is no RData file  ==>> return NULL',warning = T)
    return(NULL)
  }

  # load rdata file
  filename_variable <- load(filename)

  # check the number of loaded variables
  # if only one ==>> return
  # if multiple ==>> warning and stop
  if(length(filename_variable)==1){
    return(get(filename_variable))
  } else{
    smd_print('The RData file contains multiple variables ==>> return NULL',warning = T)
    return(NULL)
  }
}


#############################################################################
# GET FILE PATH (create required folder(s) if they are not present yet)
#############################################################################

#' @title Get (and create) folder path
#'
#' @param ... the directory name(s)
#' @param .verbose boolean to print a message if a new folder is created (default: TRUE)
#'
#' @keywords external
#' @export
smd_file_path <- function(...,.verbose=TRUE){

  # get file path based on given parameters
  file_path <- file.path(...)

  # get dirname if file_path contains an extension
  if(grepl('\\.',file_path)){
    file_dir <- dirname(file_path)
  } else{
    file_dir <- file_path
  }

  # check if the requested folder exist, and create the folder(s) if not
  if(!dir.exists(file_dir)){
    dir.create(file_dir,recursive = T)
    if(.verbose) {
      smd_print('Create directory:',file_dir)
      }
  }

  return(file_path)
}
