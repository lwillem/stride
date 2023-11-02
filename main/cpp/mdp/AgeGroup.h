/*
 *  This is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *  The software is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  You should have received a copy of the GNU General Public License
 *  along with the software. If not, see <http://www.gnu.org/licenses/>.
 *
 *  Copyright 2021
 */

/**
 * @file
 * Helpers for age groups.
 */

#pragma once

#include <vector>


namespace stride {

/// Age groups enumeration
/// TODO: refactor this to pop/Age.h?
enum AgeGroup {
    /// Age 0-4
    children,
    /// Age 5-18
    youngsters,
    /// Age 19-25
    young_adults,
    /// Age 26-64
    adults,
    /// Age 65+
    elderly
};

/// All the age groups
static const std::vector<AgeGroup> AllAgeGroups = { children, youngsters, young_adults, adults, elderly };


/// Get the age group for a given age
inline AgeGroup GetAgeGroup(unsigned int age) {
        /// Age 0-4
        if (age < 5) { return children; }
        /// Age 5-18
        else if (age < 19) { return youngsters; }
        /// Age 19-25
        else if (age < 26) { return young_adults; }
        /// Age 26-64
        else if (age < 65) { return adults; }
        /// Age 65+
        else { return elderly; }
}

enum ChildlessAgeGroup {
    /// Children
    children_c = -1,
    /// Age 18-25
    young_adults_c,
    /// Age 26-64
    adults_c,
    /// Age 65+
    elderly_c
};

/// All the age groups
static const std::vector<ChildlessAgeGroup> AllChildlessAgeGroups = { young_adults_c, adults_c, elderly_c };


/// Get the age group for a given age
inline ChildlessAgeGroup GetChildlessAgeGroup(unsigned int age) {
    /// Age 0-17 children
    if (age < 18) { return children_c; }
    /// Age 19-25
    else if (age < 26) { return young_adults_c; }
    /// Age 26-64
    else if (age < 65) { return adults_c; }
    /// Age 65+
    else { return elderly_c; }
}

} // namespace stride
