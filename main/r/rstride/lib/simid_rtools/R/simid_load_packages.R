#############################################################################
# LOAD (AND INSTALL) PACKAGES
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

#' @title Load the given packages (and install them if not present yet)
#'
#' @param ...  all packages to load
#' @keywords external
#' @export
smd_load_packages <- function(...)
{

  # collect all given package names
  all_packages <- c(...)

  # loop over the given package names
  for(i_package in all_packages){

    # if not present => install
    if(!i_package %in% rownames(installed.packages())){
      install.packages(i_package)
    }

    # load package
    library(i_package,
            character.only=TRUE,
            quietly = TRUE,
            warn.conflicts = FALSE,
            verbose = FALSE)
  }

}

