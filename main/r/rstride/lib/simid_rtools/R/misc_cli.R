#############################################################################
# COMMAND LINE INTERFACE TOOLS
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

#' @title Print (warning) message to Console
#'
#' @description Command line interface: print message with time tag, if called
#' by the master node or first parallel slave
#'
#' @param ... (parts of the) message to print
#' @param WARNING  boolean, to print the message in red
#' @param FORCED   boolean, to print the message irrespectively of the parallel thread
#'
#' @export
smd_print <- function(..., WARNING=F, FORCED=F) {

  # get function arguments
  function_arguments <- as.list(match.call(expand.dots=FALSE))$...

  # get function-call environment (to retrieve variable from that environment)
  pf <- parent.frame()

  #parse list => make character vector
  f_out <- ' '
  for(i in 1:length(function_arguments)){
    f_out <- cbind(f_out,eval(unlist(function_arguments[[i]]),envir = pf))
  }

  # add a space to each function arguments
  function_arguments <- paste(f_out,collapse = ' ')

  # set text color: default (black/white) or red (warning)
  text_color         <- smd_get_console_color(WARNING)

  # print time + arguments (without spaces)
  cli_out <- paste0(c('echo "',text_color, '[',format(Sys.time(),'%H:%M:%S'),']',
                      function_arguments,'"'),collapse = '')

  # print if function is called by master-node or first slave
  if(!exists('par_nodes_info') ||
     Sys.getpid() == par_nodes_info$pid_master ||
     FORCED){
      system(cli_out)
  }

  # add to R warnings
  if(WARNING){
    cli_warning <- paste0(c(text_color, '[',format(Sys.time(),'%H:%M:%S'),']',
                            function_arguments),collapse = '')
    warning(cli_warning,
            call. = FALSE, immediate.=FALSE)
  }
}


#' @title Print progress bar
#'
#' @description Progress bar interface
#'
#' @param i_current       current iteration
#' @param i_total         total iterations
#' @param time_stamp_loop starting time
#' @param par_nodes_info  info about parallel nodes
#'
#' @export
smd_print_progress <- function(i_current,i_total, time_stamp_loop = Sys.time(),par_nodes_info = NA){

  # print if function is called by first node
  if(!is.na('par_nodes_info') && (Sys.getpid() == par_nodes_info$pid_master || Sys.getpid() == par_nodes_info$pid_slave1)){

    # calculate progress
    progress_scen        <- floor(i_current/i_total*100)
    progress_time        <- round(difftime(Sys.time(),time_stamp_loop,units = "min"),digits=1)

    # estimate remaining time (after 15%)
    time_label <- ''
    if(progress_time > 0 & progress_scen > 15 & progress_scen < 99) {
      estim_time <-  round(progress_time / progress_scen * (100-progress_scen),digits=1)
      if(estim_time<1) {estim_time <- '<1'}
      time_label <- paste0('[',estim_time,' min remaining]')
     }

    smd_print('RUNNING...',i_current,'/',i_total,time_label,FORCED=TRUE)
  }
}

#' @title Print progress bar [DEPRECATED]
#'
#' @description please use smd_print_progress()
#' @export
smd_progress <- function(i_current,i_total, time_stamp_loop = Sys.time(),par_nodes_info = NA){

  smd_print('PLEASE USE smd_print_progress() insead of smd_progress()',WARNING = T)

}

#' @title Get default Rstudio console color
#'
#' @description The default color of the console (or terminal) depends on whether R(studio)
#' is used in DARK mode. This function return 'black' or 'white', based on the current theme.
#' If the function is used in R, an empty color is returned.
#'
#' @keywords internal
smd_get_console_color <- function(WARNING=FALSE){

  if(!rstudioapi::isAvailable()){
    return('')
  }

  if(WARNING){
    return('\033[0;31m') #red
  }

  if(rstudioapi::getThemeInfo()$dark){
    return('\033[0;37m') # white
  } else{
    return('\033[0;30m') # black
  }

}
