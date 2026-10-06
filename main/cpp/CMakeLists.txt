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


#============================================================================
# Sources.
#============================================================================
set(STRIDE_SRC
    #---
    calendar/Calendar.cpp
    #---
    contact/AgeContactProfile.cpp
    contact/ContactDivider.cpp
    contact/ContactHeterogeneitySeeder.cpp
    contact/ContactPool.cpp
    contact/ContactPoolSys.cpp
    contact/ContactProbabilityRule.cpp
    contact/ContactType.cpp
    contact/EventLogMode.cpp
    contact/Infector.cpp
    contact/PoolCharacteristicsSeeder.cpp
    #--
    execs/SimController.cpp
    #---
    health/DiseaseSeeder.cpp
    health/Health.cpp
    health/HealthSeeder.cpp
    health/ImmunitySeeder.cpp
    health/TransmissionProfile.cpp
    #---
    healthcare/PublicHealthAgency.cpp
    #---
    pop/Person.cpp
    pop/Population.cpp
    pop/PopBuilder.cpp
    pop/PopFileLayout.cpp
    pop/SurveyManager.cpp
    #---
    sim/Sim.cpp
    sim/SimBuilder.cpp
    #---
    util/ConfigInfo.cpp
    util/Date.cpp
    util/FileSys.cpp
    util/GammaDistribution.cpp
    util/LogUtils.cpp
    util/PopSnapshotWriter.cpp
    util/Ptree.cpp
    util/RunConfigManager.cpp
    util/SummaryFile.cpp
    util/TimeStamp.cpp
    util/TimeToString.cpp
)

#============================================================================
# Build the Stride library.
#============================================================================
add_library(libstride ${STRIDE_SRC})
target_compile_definitions(libstride PUBLIC "-DSTRIDE_GIT_HASH=${STRIDE_WC_REVISION_HASH}")
target_compile_definitions(libstride PUBLIC "-DPROCCOUNT=${PROCCOUNT}")
target_link_libraries(libstride PRIVATE ${LIBS})

# CMakeCPP.cmake adds ${OpenMP_CXX_FLAGS} to the COMPILE flags but never linked the
# runtime, so a successful detection would still fail at link time. The imported target
# carries both halves. Guarded, so this is inert whenever OpenMP is not found and the
# dummy stubs in main/resources/lib/domp are used instead.
if(OPENMP_FOUND AND TARGET OpenMP::OpenMP_CXX)
    target_link_libraries(libstride PUBLIC OpenMP::OpenMP_CXX)
endif()
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
