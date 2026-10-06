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

#============================================================================
# UNZIP files from data/* to ${DATA_INSTALL_LOCATION}
#============================================================================
# Create build folder to unzip all compressed files
set(UNZIP_FOLDER ${CMAKE_BINARY_DIR}/main/resources/data_zip.noindex/)
execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${UNZIP_FOLDER})

# Unzip all ZIP files.
# CONFIGURE_DEPENDS so that an archive added to data/ after this build directory was
# created still gets unpacked: without it the glob is evaluated once, the unzip below
# never re-runs, and the archive silently never reaches the install (plan F7.1).
file(GLOB ZIP_FILES CONFIGURE_DEPENDS "data/*.zip")
# Unpack an archive only when it is newer than its stamp (or has none yet), so a
# reconfigure does not rewrite ~650MB of extracted files each time. The stamps live
# outside UNZIP_FOLDER so the install glob below does not pick them up. Delete a
# stamp (or the build directory) to force a fresh unpack.
set(UNZIP_STAMP_FOLDER ${CMAKE_BINARY_DIR}/main/resources/data_zip_stamps.noindex/)
execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${UNZIP_STAMP_FOLDER})
foreach( file_i ${ZIP_FILES} )
	get_filename_component(name_i ${file_i} NAME)
	set(stamp_i ${UNZIP_STAMP_FOLDER}${name_i}.stamp)
	if(${file_i} IS_NEWER_THAN ${stamp_i})
		execute_process(COMMAND ${CMAKE_COMMAND} -E tar xzf ${file_i}
	                    WORKING_DIRECTORY ${UNZIP_FOLDER}
	                    RESULT_VARIABLE unzip_result_i)
		if(unzip_result_i EQUAL 0)
			file(TOUCH ${stamp_i})
		else()
			message(WARNING "Failed to unpack ${file_i}")
		endif()
	endif()
endforeach()
unset(UNZIP_STAMP_FOLDER)

# Install unzipped files.
# No CONFIGURE_DEPENDS here: this globs the build directory, which the loop above keeps
# in step with data/*.zip in this same configure pass, so it is current by construction.
file(GLOB DEPLOY_FILES ${UNZIP_FOLDER}*)
if(NOT STRIDE_INSTALL_POPULATION_PDF)
	file(GLOB DEPLOY_PDF_FILES ${UNZIP_FOLDER}*.pdf)
	list(REMOVE_ITEM DEPLOY_FILES ${DEPLOY_PDF_FILES})
	unset(DEPLOY_PDF_FILES)
endif()
INSTALL(FILES ${DEPLOY_FILES} DESTINATION ${DATA_INSTALL_LOCATION})
unset(UNZIP_FOLDER)

#============================================================================
# Copy all non-zipped files from data/* to ${DATA_INSTALL_LOCATION}
#============================================================================
file(GLOB DEPLOY_FILES CONFIGURE_DEPENDS "data/*")
list(REMOVE_ITEM DEPLOY_FILES ${ZIP_FILES})
if(NOT STRIDE_INSTALL_POPULATION_PDF)
	file(GLOB DEPLOY_PDF_FILES CONFIGURE_DEPENDS "data/*.pdf")
	list(REMOVE_ITEM DEPLOY_FILES ${DEPLOY_PDF_FILES})
	unset(DEPLOY_PDF_FILES)
endif()
INSTALL(FILES ${DEPLOY_FILES} DESTINATION ${DATA_INSTALL_LOCATION})
unset(ZIP_FILES)
unset(DEPLOY_FILES)

#============================================================================
# Copy all files from config/* to ${CONFIG_INSTALL_LOCATION}
#============================================================================
file(GLOB DEPLOY_FILES CONFIGURE_DEPENDS "config/*")
INSTALL(FILES ${DEPLOY_FILES} DESTINATION ${CONFIG_INSTALL_LOCATION})
unset(DEPLOY_FILES)

#============================================================================
# Copy all files from rstride_test/* to ${TESTS_INSTALL_LOCATION}
#============================================================================
file(GLOB DEPLOY_FILES CONFIGURE_DEPENDS "rstride_test/*")
INSTALL(FILES ${DEPLOY_FILES} DESTINATION ${TESTS_INSTALL_LOCATION})
unset(DEPLOY_FILES)


#============================================================================
# Add subdirs:
#============================================================================
add_subdirectory(lib)

#############################################################################