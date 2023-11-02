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
#  Copyright 2021
#############################################################################


#============================================================================
# Sources.
#============================================================================
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(pybind11.cmake)
# Stride MDP code
set(PYLIBSTRIDE_SRC
    "../../cpp/mdp/MDP.cpp"
    "../../cpp/mdp/MDPRunner.cpp"
    )

# Pybind wrappers
set(PYBIND_SRC
    "mdp/MDP.cpp"
    "mdp/AgeGroup.cpp"
    "mdp/Vaccines.cpp"
#    "mdp/Health.cpp"
    "libstride.cpp"
    )

pybind11_add_module(pylibstride ${PYLIBSTRIDE_SRC} ${PYBIND_SRC})

target_compile_definitions(pylibstride PUBLIC "-DSTRIDE_GIT_HASH=${STRIDE_WC_REVISION_HASH}")
target_compile_definitions(pylibstride PUBLIC "-DPROCCOUNT=${PROCCOUNT}")
target_link_libraries(pylibstride PRIVATE libstride)
set_target_properties(pylibstride
        PROPERTIES
        POSITION_INDEPENDENT_CODE ON
        PREFIX ""
        )

#============================================================================
# Build the Python Stride library.
#============================================================================
install(TARGETS pylibstride COMPONENT python LIBRARY DESTINATION "${PYTHON_LIBRARY_DIR}")
