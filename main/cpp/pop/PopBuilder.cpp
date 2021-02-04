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
 *  Copyright 2017, 2018, Kuylen E, Willem L, Broeckhove J
 *  Copyright 2018, 2019 Jan Broeckhove and Bistromatics group.
 */

/**
 * @file
 * Initialize populations: implementation.
 */

#include "contact/ContactType.h"
#include "contact/IdSubscriptArray.h"
#include "pop/Population.h"
#include "pop/SurveySeeder.h"
#include "util/FileSys.h"
#include "util/RnMan.h"
#include "util/StringUtils.h"
#include "util/LogUtils.h"

#include <boost/property_tree/ptree.hpp>
#include <fstream>
#include "PopBuilder.h"
#include "../contact/EventLogMode.h"

namespace stride {

using namespace ContactType;

using namespace util;
using namespace boost::property_tree;
using namespace std;

PopBuilder::PopBuilder(const boost::property_tree::ptree& config,
                                       std::shared_ptr<spdlog::logger> strideLogger)
    : m_config(config), m_stride_logger(std::move(strideLogger))
{
        if (!m_stride_logger) {
                m_stride_logger = util::LogUtils::CreateNullLogger("PopBuilder_logger");
        }
}

shared_ptr<Population> PopBuilder::MakePersons(shared_ptr<Population> pop)
{
        //------------------------------------------------
        // Read persons from file.
        //------------------------------------------------
        const auto fileName = m_config.get<string>("run.population_file");
        m_stride_logger->info("Building default population from file {}.", fileName);

        const auto use_install_dirs = m_config.get<bool>("run.use_install_dirs");
        const auto filePath         = (use_install_dirs) ? FileSys::GetDataDir() /= fileName : filesys::path(fileName);
        if (!is_regular_file(filePath)) {
                throw runtime_error(string(__func__) + "> Population file " + filePath.string() + " not present.");
        }

        ifstream popFile;
        popFile.open(filePath.string());
        if (!popFile.is_open()) {
                throw runtime_error(string(__func__) + "> Error opening population file " + filePath.string());
        }

        // get age break between 2 school types
        //TODO: rename school types and/or add 3rd for secondary school
        const unsigned int age_break_school_types = m_config.get<unsigned int>("run.age_break_school_types",18);

        string line;
        getline(popFile, line); // step over file header
        auto headers   = Split(line, ",");

        unsigned int person_id = 0U;

        while (getline(popFile, line)) {
                const auto values               = Split(line, ",");
                const auto age                  = FromString<unsigned int>(values[0]);
                const auto householdId          = FromString<unsigned int>(values[1]);
                auto schoolId                   = FromString<unsigned int>(values[2]);
                const auto workId               = FromString<unsigned int>(values[3]);
                const auto primaryCommunityId   = FromString<unsigned int>(values[4]);
                const auto secondaryCommunityId = FromString<unsigned int>(values[5]);

                unsigned int householdClusterId = 0;
                if(values.size() == 7 && Trim(ToString(headers[6]),ToString('"')) == "household_cluster_id"){
                	householdClusterId = FromString<unsigned int>(values[6]);
                }

                unsigned int collectivityId = 0;

				if(values.size() == 7 && Trim(ToString(headers[6]),ToString('"')) == "collectivity_id"){
					collectivityId = FromString<unsigned int>(values[6]);
				}

                //TODO: rename school types to current approach
                unsigned int collegeId = 0;
                if(schoolId != 0 && age >= age_break_school_types && age < 23){
                	collegeId = schoolId;
                	schoolId = 0;
                }

                pop->CreatePerson(person_id, age, householdId, schoolId, collegeId, workId, primaryCommunityId,
                                  secondaryCommunityId, householdClusterId, collectivityId);
                ++person_id;
        }

        popFile.close();

        m_stride_logger->trace("Done building default population.");

        return pop;
}

shared_ptr<Population> PopBuilder::Build(shared_ptr<Population> pop)
{
        //------------------------------------------------
        // Add persons
        //------------------------------------------------
        MakePersons(pop);

        // --------------------------------------------------------------
        // Determine maximum pool ids in population.
        // --------------------------------------------------------------
        IdSubscriptArray<unsigned int> maxIds{0U};
        for (const auto& p : *pop) {
                for (Id typ : IdList) {
                        maxIds[typ] = max(maxIds[typ], p.GetPoolId(typ));
                }
        }
        // --------------------------------------------------------------
        // Initialize poolSys with empty ContactPools (even for Id=0).
        // --------------------------------------------------------------
        for (Id typ : IdList) {
                for (unsigned int i = 1; i < maxIds[typ] + 1; i++) {
                        pop->RefPoolSys().CreateContactPool(typ);
                }
        }

        // --------------------------------------------------------------
        // Insert persons (pointers) in their contactpools. Having Id 0
        // means "not belonging pool of that type" (e.g. school/ work -
        // cannot belong to both, or e.g. out-of-work).
        //
        // Pools are uniquely identified by (type, subscript) and a Person
        // belongs, per type, to the pool with subscript p.GetPoolId(type).
        // Defensive measure: we have a pool for Id 0 and leave it empty.
        // --------------------------------------------------------------
        for (auto& p : *pop) {
                for (Id typ : IdList) {
                        const auto poolId = p.GetPoolId(typ);
                        if (poolId > 0) {
                                pop->RefPoolSys().RefPools(typ)[poolId].AddMember(&p);
                        }
                }
        }


        return pop;
}

} // namespace stride
