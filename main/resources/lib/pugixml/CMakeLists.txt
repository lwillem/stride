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
#  Copyright 2026
#############################################################################
#
#  Vendored pugixml (MIT licensed), used as a lightweight replacement for
#  boost::property_tree's XML parsing/writing.
#
#############################################################################

add_library(pugixml src/pugixml.cpp)
target_include_directories(pugixml PUBLIC src)
set_target_properties(pugixml
    PROPERTIES
        POSITION_INDEPENDENT_CODE ON
)
