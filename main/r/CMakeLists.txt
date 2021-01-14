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
#  Copyright 2020, Willem L, Kuylen E & Broeckhove J
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
	PATTERN "*.zip" EXCLUDE
)


# rStride analysis scripts
INSTALL( FILES 			
            rStride_abc.R
            rStride_all.R
            rStride_contacts.R
            rStride_create_config.R
            rStride_ensemble.R
            rStride_param.R
            rStride_parse_output.R
            rStride_explore.R
            rStride_gtester_covid19.R
            rStride_r0.R 	
            rStride_intervention_baseline.R
            rStride_intervention.R
            rStride_superspreading.R	
            rStride_trans_prob_estimation.R			 	 
            rStride_universal_testing.R
            rStride_pop.R				 

	DESTINATION ${BIN_INSTALL_LOCATION}  
	PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_WRITE GROUP_READ
)

#============================================================================
# UNZIP R-packages from rstride/lib to ${BIN_INSTALL_LOCATION}/bin/rstride/lib
#============================================================================
# Create build folder to unzip all compressed files
set(UNZIP_FOLDER ${CMAKE_BINARY_DIR}/main/r/package_zip/)
execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${UNZIP_FOLDER})

# Unzip all ZIP files
file(GLOB ZIP_FILES "rstride/lib/*.zip")
foreach( file_i ${ZIP_FILES} )
	execute_process(COMMAND ${CMAKE_COMMAND} -E tar xzf ${file_i}
                    WORKING_DIRECTORY ${UNZIP_FOLDER})
endforeach()

# Install unzipped files

INSTALL(DIRECTORY ${UNZIP_FOLDER}
	    DESTINATION ${BIN_INSTALL_LOCATION}/rstride/lib
)
unset(UNZIP_FOLDER)

#############################################################################
