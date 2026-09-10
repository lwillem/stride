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
 * Implementation for the PopSnapshotWriter class.
 */

#include "util/PopSnapshotWriter.h"

#include "contact/ContactType.h"
#include "pop/Age.h"
#include "pop/Person.h"
#include "util/FileSys.h"

#include <fstream>
#include <map>
#include <sstream>
#include <vector>

namespace stride {
namespace util {

void PopSnapshotWriter::Write(const stride::util::ptree& config, std::shared_ptr<Population> pop)
{
        const auto outputPrefix = config.get<std::string>("run.output_prefix");
        WriteHouseholdMemberships(outputPrefix, pop);
        WriteSusceptiblesByAge(outputPrefix, pop);
}

void PopSnapshotWriter::WriteHouseholdMemberships(const std::string& outputPrefix, std::shared_ptr<Population> pop)
{
        struct HouseholdData
        {
                std::vector<int> ages;
                std::vector<int> susceptible;
        };
        std::map<unsigned int, HouseholdData> households;

        for (const auto& p : *pop) {
                const auto householdId = p.GetPoolId(ContactType::Id::Household);
                auto&      entry       = households[householdId];
                entry.ages.push_back(static_cast<int>(p.GetAge()));
                entry.susceptible.push_back(p.GetHealth().IsSusceptible() ? 1 : 0);
        }

        const auto    path = FileSys::BuildPath(outputPrefix, "households.csv");
        std::ofstream csvFile(path.string());
        csvFile << "household_id,ages,susceptible\n";

        for (const auto& kv : households) {
                std::ostringstream agesStr;
                agesStr << "[";
                for (size_t i = 0; i < kv.second.ages.size(); ++i) {
                        agesStr << kv.second.ages[i];
                        if (i + 1 < kv.second.ages.size()) {
                                agesStr << ", ";
                        }
                }
                agesStr << "]";

                std::ostringstream suscStr;
                suscStr << "[";
                for (size_t i = 0; i < kv.second.susceptible.size(); ++i) {
                        suscStr << kv.second.susceptible[i];
                        if (i + 1 < kv.second.susceptible.size()) {
                                suscStr << ", ";
                        }
                }
                suscStr << "]";

                csvFile << kv.first << ",\"" << agesStr.str() << "\",\"" << suscStr.str() << "\"\n";
        }
}

void PopSnapshotWriter::WriteSusceptiblesByAge(const std::string& outputPrefix, std::shared_ptr<Population> pop)
{
        const unsigned int       maxAge = MaximumAge();
        std::vector<unsigned int> susceptibleCount(maxAge + 1, 0U);
        std::vector<unsigned int> immuneCount(maxAge + 1, 0U);

        for (const auto& p : *pop) {
                const auto age          = static_cast<unsigned int>(p.GetAge());
                const auto effectiveAge = (age <= maxAge) ? age : maxAge;
                if (p.GetHealth().IsSusceptible()) {
                        susceptibleCount[effectiveAge]++;
                } else {
                        immuneCount[effectiveAge]++;
                }
        }

        const auto    path = FileSys::BuildPath(outputPrefix, "susceptibles_by_age.csv");
        std::ofstream csvFile(path.string());
        csvFile << "age,susceptible,immune\n";
        for (unsigned int age = 0; age <= maxAge; ++age) {
                csvFile << age << "," << susceptibleCount[age] << "," << immuneCount[age] << "\n";
        }
}

} // namespace util
} // namespace stride
