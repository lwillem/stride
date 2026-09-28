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
#  Copyright 2024
#############################################################################

#================================================================================================
# Copy the rStride scripts to the bin folder.
#================================================================================================
# rStride main files
INSTALL(DIRECTORY
	   rstride
	DESTINATION ${BIN_INSTALL_LOCATION}
	PATTERN "*.Rproj*" EXCLUDE
	PATTERN "*.Rhistory" EXCLUDE
)

# rStride analysis scripts
#
# Globbed rather than enumerated by hand. The list used to be maintained manually and
# had drifted: rStride_measles_default_param.R, rStride_param.R and rStride_r0_measles.R
# existed in the repository but never reached bin/, so scripts sourcing them failed with
# a missing-file error that pointed at the install directory rather than at this list.
#
# CONFIGURE_DEPENDS makes CMake re-evaluate the glob at build time, so adding a script is
# just adding the file. Note the top-level Makefile runs `configure` before `all`, so a
# plain `make install` picks up new scripts either way.
FILE(GLOB RSTRIDE_SCRIPTS CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/rStride_*.R")
LIST(SORT RSTRIDE_SCRIPTS)

INSTALL( FILES ${RSTRIDE_SCRIPTS}
	DESTINATION ${BIN_INSTALL_LOCATION}  
	PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_WRITE GROUP_READ
)

#############################################################################
