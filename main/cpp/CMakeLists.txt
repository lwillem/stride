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
#  Copyright 2017, 2018 Willem L, Kuylen E & Broeckhove J
#############################################################################


#============================================================================
# Sources.
#============================================================================
set(STRIDE_SRC
    #---
    calendar/Calendar.cpp
    #---
    contact/AgeContactProfile.cpp
    contact/ContactPool.cpp
    contact/ContactPoolSys.cpp
    contact/ContactType.cpp
    contact/EventLogMode.cpp
    contact/Infector.cpp
    contact/NonComplianceSeeder.cpp
    #---
    disease/DiseaseSeeder.cpp
    disease/Health.cpp
    disease/HealthSeeder.cpp
    disease/ImmunitySeeder.cpp
    disease/PublicHealthAgency.cpp
    disease/TransmissionProfile.cpp
    disease/UniversalTesting.cpp
    #---
    execs/ControlHelper.cpp
    execs/SimController.cpp
    #---
    pop/Person.cpp
    pop/Population.cpp
    pop/PopBuilder.cpp
    pop/SurveySeeder.cpp
    #---
    sim/SimRunner.cpp
    sim/Sim.cpp
    sim/SimBuilder.cpp
    sim/HospitalisationConfig.cpp
    sim/event/Id.cpp
    sim/event/Subject.cpp
    #---
    util/ConfigInfo.cpp
    util/CSV.cpp
    util/CSVRow.cpp
    util/FileSys.cpp
    util/LogUtils.cpp
    util/Rn.cpp
    util/RnMan.cpp
    util/RunConfigManager.cpp
    util/TimeStamp.cpp
    util/TimeToString.cpp
    #---
    viewers/CliViewer.cpp
    viewers/InfectedFile.cpp
    viewers/InfectedFileViewer.cpp
    viewers/SummaryFile.cpp
    viewers/SummaryFileViewer.cpp
)

#============================================================================
# Build the Stride library.
#============================================================================
add_library(libstride ${STRIDE_SRC})
target_compile_definitions(libstride PUBLIC "-DSTRIDE_GIT_HASH=${STRIDE_WC_REVISION_HASH}")
target_compile_definitions(libstride PUBLIC "-DPROCCOUNT=${PROCCOUNT}")
target_link_libraries(libstride PRIVATE ${LIBS})
set_target_properties(libstride
    PROPERTIES
        POSITION_INDEPENDENT_CODE ON
        PREFIX ""
)
include(CheckIPOSupported)

# For future use: inter procedural optimization
# Disabled for now because of issues with gcc 8.1.0 on travis.
#check_ipo_supported(RESULT ipo_result)
#if(ipo_result)
#    set_target_properties(libstride PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
#endif()
#message( STATUS "CheckIPOSupported: ${ipo_result}" )

#============================================================================
# Build the stride binary.
#============================================================================
add_executable(stride execs/stride.cpp)
target_link_libraries(stride libstride)
install(TARGETS stride DESTINATION ${BIN_INSTALL_LOCATION})

#############################################################################
