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
 *  Copyright 2026
 */

/**
 * @file
 * Header for the PopSnapshotWriter class.
 */

#pragma once

#include "pop/Population.h"
#include "util/Ptree.h"

#include <memory>
#include <string>

namespace stride {
namespace util {

/**
 * Writes a population-wide snapshot (household memberships & susceptibles-by-age)
 * to CSV files. Intended to be called once, at the start of a simulation, to record
 * who is initially infected/susceptible/immunized.
 */
class PopSnapshotWriter
{
public:
        /// Write both households.csv and susceptibles_by_age.csv to the run's output prefix.
        static void Write(const stride::util::ptree& config, std::shared_ptr<Population> pop);

private:
        /// Write households.csv: one row per household, listing member ages & susceptibility.
        static void WriteHouseholdMemberships(const std::string& outputPrefix, std::shared_ptr<Population> pop);

        /// Write susceptibles_by_age.csv: susceptible/immune counts per age.
        static void WriteSusceptiblesByAge(const std::string& outputPrefix, std::shared_ptr<Population> pop);
};

} // namespace util
} // namespace stride
