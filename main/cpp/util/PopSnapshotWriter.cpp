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
#include "util/StringUtils.h"

#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace stride {
namespace util {

void PopSnapshotWriter::Write(const stride::util::ptree& config, std::shared_ptr<Population> pop)
{
        const auto outputPrefix = config.get<std::string>("run.output_prefix");
        WriteHouseholdMemberships(outputPrefix, pop);
        WriteSusceptiblesByAge(outputPrefix, pop);
        WritePopulationSnapshot(config, outputPrefix, pop);
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
                // "susceptible" here means "not immune" (natural or vaccine-induced),
                // per Person::IsImmune() -- see class comment in the header.
                entry.susceptible.push_back(p.IsImmune() ? 0 : 1);
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
        const unsigned int        maxAge = MaximumAge();
        std::vector<unsigned int> susceptibleCount(maxAge + 1, 0U);
        std::vector<unsigned int> immuneCount(maxAge + 1, 0U);

        for (const auto& p : *pop) {
                const auto age          = static_cast<unsigned int>(p.GetAge());
                const auto effectiveAge = (age <= maxAge) ? age : maxAge;
                // "susceptible" here means "not immune" (natural or vaccine-induced),
                // per Person::IsImmune() -- see class comment in the header.
                if (p.IsImmune()) {
                        immuneCount[effectiveAge]++;
                } else {
                        susceptibleCount[effectiveAge]++;
                }
        }

        const auto    path = FileSys::BuildPath(outputPrefix, "susceptibles_by_age.csv");
        std::ofstream csvFile(path.string());
        csvFile << "age,susceptible,immune\n";
        for (unsigned int age = 0; age <= maxAge; ++age) {
                csvFile << age << "," << susceptibleCount[age] << "," << immuneCount[age] << "\n";
        }
}

void PopSnapshotWriter::WritePopulationSnapshot(const stride::util::ptree& config, const std::string& outputPrefix,
                                                std::shared_ptr<Population> pop)
{
        // --------------------------------------------------------------
        // Re-read just the header line of the population input file, using the
        // same parsing rules as PopBuilder::MakePersons, so our output mirrors
        // whatever column layout that file actually used (with/without a
        // profession column, with/without household_cluster_id / collectivity_id).
        // --------------------------------------------------------------
        const auto    fileName = config.get<std::string>("run.population_file");
        std::ifstream popFile(fileName);
        if (!popFile.is_open()) {
                throw std::runtime_error(std::string(__func__) + "> Error opening population file " + fileName);
        }
        std::string headerLine;
        std::getline(popFile, headerLine);
        popFile.close();

        const bool        useSemicolon = headerLine.find(';') != std::string::npos;
        const std::string sep          = useSemicolon ? ";" : ",";
        const auto        headers      = Split(headerLine, sep);

        const bool         hasProfession = headers.size() > 2 && Trim(ToString(headers[2]), ToString('"')) == "worker";
        const unsigned int professionAdj = hasProfession ? 2 : 0;
        const bool         hasExtraColumn = headers.size() == (7 + professionAdj);

        std::string extraId;
        if (hasExtraColumn) {
                extraId = Trim(ToString(headers[6 + professionAdj]), ToString('"'));
        }
        const bool hasHouseholdClusterId = extraId == "household_cluster_id";
        const bool hasCollectivityId     = extraId == "collectivity_id";

        // --------------------------------------------------------------
        // Write population_snapshot.csv: original columns + immunity_status.
        // --------------------------------------------------------------
        const auto    path = FileSys::BuildPath(outputPrefix, "population_snapshot.csv");
        std::ofstream csvFile(path.string());

        csvFile << "age" << sep;
        if (hasProfession) {
                csvFile << "person_id" << sep << "profession" << sep;
        }
        csvFile << "household_id" << sep << "school_id" << sep << "workplace_id" << sep << "community_weekend_id"
                 << sep << "community_weekday_id" << sep;
        if (hasHouseholdClusterId) {
                csvFile << "household_cluster_id" << sep;
        } else if (hasCollectivityId) {
                csvFile << "collectivity_id" << sep;
        }
        csvFile << "immunity_status\n";

        for (const auto& p : *pop) {
                csvFile << static_cast<int>(p.GetAge()) << sep;
                if (hasProfession) {
                        csvFile << p.GetId() << sep << p.GetProfession() << sep;
                }
                csvFile << p.GetPoolId(ContactType::Id::Household) << sep
                         << p.GetPoolId(ContactType::Id::School) << sep
                         << p.GetPoolId(ContactType::Id::Workplace) << sep
                         << p.GetPoolId(ContactType::Id::CommunityWeekend) << sep
                         << p.GetPoolId(ContactType::Id::CommunityWeekday) << sep;
                if (hasHouseholdClusterId) {
                        csvFile << p.GetPoolId(ContactType::Id::HouseholdCluster) << sep;
                } else if (hasCollectivityId) {
                        csvFile << p.GetPoolId(ContactType::Id::Collectivity) << sep;
                }
                // "immune" means Person::IsImmune() -- natural or vaccine-induced immunity
                // assigned by ImmunitySeeder -- not Health's disease-progression status.
                csvFile << (p.IsImmune() ? "immune" : "susceptible") << "\n";
        }
}

} // namespace util
} // namespace stride