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
 *  Copyright 2026, Kuylen E, Willem L, Broeckhove J
 */

/**
 * @file
 * Header for the ContactProbabilityRule.
 */

#pragma once

#include <string>

namespace stride {

/**
 * Enum specifying how the two age-specific contact probabilities of a candidate
 * contact pair are combined into a single contact probability:
 * \li Min:  the minimum of both probabilities (default)
 * \li Mean: the average of both probabilities
 *
 * NOTE: this rule is a calibration dependency. The <transmission> b0/b1/b2 fit in
 * every disease configuration file maps a requested R0 onto a transmission
 * probability, and that mapping is only valid for the rule it was fitted under.
 * Changing the rule without refitting silently changes the realised R0.
 */
namespace ContactProbabilityRule {

enum class Id
{
        Min  = 0U,
        Mean = 1U
};

/// Converts a ContactProbabilityRule value to the corresponding name.
std::string ToString(ContactProbabilityRule::Id r);

/// Check whether the string is the name of a ContactProbabilityRule value.
bool IsRule(const std::string& s);

/// Converts a string with a name to a ContactProbabilityRule value.
ContactProbabilityRule::Id ToRule(const std::string& s);

} // namespace ContactProbabilityRule
} // namespace stride
