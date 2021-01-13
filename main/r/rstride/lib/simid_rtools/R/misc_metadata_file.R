#############################################################################
# METADATA FILE CREATION AND MANIPULATION
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


#' @title Create a text file with METADATA
#'
#' @description Create a file with METADATA on e.g. current software state, run info, and other info.
#'
#' @param output_dir  the directory where the new file will be created.
#' @param file_prefix a prefix for the METADATA file (optional, but overruled by 'data_file').
#' @param run_tag     the tag to specify the model run (optional).
#' @param run_time    the run time (optional).
#' @param root_dir    the root directory, with a DESCRIPTION or METADATA file to start from (default = '.')
#' @param other_info  additional info to include (optional).
#' @param data_file   the corresponding DATA file, to create a correspoding file name for the METADATA file (optional).
#' @param main_title  new 'Title' of the METADATA file (optional).
#' @param sub_title   subtitle for the new section to add, which will be capitalized (default: PROCESSING DETAILS).
#'
#'
#' @export
smd_create_metadata_file <- function(output_dir, data_file = NA, run_tag = NA,
                                     run_time = NA, root_dir = '.', other_info = NA,
                                     file_prefix = '', main_title=NA, sub_title='PROCESSING DETAILS'){

  # set default metadata filename
  metadata_filename <- "METADATA.txt"

  # if data_file provided, create new name based on the data_file
  # else, if file_prexis provided, add to name
  if(!is.na(data_file)){
    # get data file name without extension
    data_file_name <- file_path_sans_ext(data_file)
    # add base name to metadata file name
    metadata_filename <-  paste(data_file_name,metadata_filename,sep='_')
    } else if(nchar(file_prefix)>0){
    metadata_filename <-  paste(file_prefix,metadata_filename,sep='_')
  }

  # set full metadata filename with output directory
  metadata_filename    <- smd_file_path(output_dir,metadata_filename)

  # look for a METADATA file in the root folder
  source_filename      <- dir(root_dir,pattern = "METADATA",full.names=T)

  # else, look for a DESCRIPTION file in the root folder
  if(length(source_filename)==0){
    source_filename      <- dir(root_dir,pattern = "DESCRIPTION",full.names=T)
  }

  # if present, load source file with version number, software info and commit id
  # else, start from scratch
  if(length(source_filename)==1){
    metadata <- readLines(source_filename)
  } else{
    metadata <- 'Title: METADATA'
  }

  if(!is.na(main_title)){
    bool_title_line <- grepl('Title:',metadata)
    metadata[bool_title_line] <- paste('Title:',main_title)
  }

  # add blank line, dashed line and subtitle (optional)
  metadata <- c(metadata,
                '',
                '----------------------------------------',
                toupper(sub_title),
                '')

  # if provided, add run tag
  if(!is.na(run_tag)){
    metadata <- c(metadata,paste('Run_tag:\t',run_tag))
  }

  # if provided, add run time
  if(!is.na(run_time)){
    metadata <- c(metadata,paste('Run_time:\t',run_time))
  }

  # add current date and some environment variables
  metadata <- c(metadata,
                paste('Run_date:\t',format(Sys.time())),
                paste('Run_R_version:\t',sessionInfo()$R.version$version.string),
                paste('Run_platform:\t',sessionInfo()$platform),
                paste('Run_workdir:\t',getwd())
  )

  # if provided, add other info
  if(any(!is.na(other_info))){
    metadata <- c(metadata,'',other_info)
  }

  # create the METADATA file in the output directory
  writeLines(metadata, metadata_filename)
}


