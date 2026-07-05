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
 *  Copyright 2019, Willem L, Kuylen E, Broeckhove J
 */

/**
 * @file
 * Contact profile.
 */

#pragma once

#include "ContactType.h"
#include "pop/Age.h"

#include "util/Ptree.h"
#include <array>

namespace stride {

/**
 * Contact rate as a function of age.
 */
class AgeContactProfile : public std::array<double, MaximumAge() + 1>
{
public:
        /// Need to keep the default constructor available.
        AgeContactProfile() = default;

        /// Explicitly initialize
        AgeContactProfile(ContactType::Id poolType, const stride::util::ptree& contactPt);
};

} // namespace stride
