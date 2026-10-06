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
#include "pop/Person.h"
#include "pop/Population.h"
#include "pop/PopFileLayout.h"
#include "pop/SurveyManager.h"
#include "util/FileSys.h"
#include "util/RnMan.h"
#include "util/StringUtils.h"
#include "util/LogUtils.h"

#include "util/Ptree.h"
#include <fstream>
#include <map>
#include "PopBuilder.h"
#include "../contact/EventLogMode.h"

using namespace stride::ContactType;
using namespace stride::util;
using namespace std;

namespace stride {

PopBuilder::PopBuilder(const stride::util::ptree& config,
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

    const filesys::path filePath{fileName};
    if (!is_regular_file(filePath)) {
        throw runtime_error(string(__func__) + "> Population file " + filePath.string() + " not present.");
    }

    ifstream popFile;
    popFile.open(filePath.string());
    if (!popFile.is_open()) {
        throw runtime_error(string(__func__) + "> Error opening population file " + filePath.string());
    }

    // The first line is either a header or, in a file without one, the first person.
    string line;
    getline(popFile, line);
    const auto layout = PopFileLayout::FromFirstLine(line);
    if (layout.HasHeader()) {
        m_stride_logger->info("Population file columns resolved by header name.");
    } else {
        m_stride_logger->warn("Population file has no header: columns read by position.");
    }
    for (const auto& c : layout.PositionalColumns()) {
        m_stride_logger->warn("Population file: unrecognised column name {}.", c);
    }
    for (const auto& c : layout.IgnoredColumns()) {
        m_stride_logger->info("Population file: column {} not used.", c);
    }

    // Read persons from file
    unsigned int default_person_id = 0U;
    const auto   add_person        = [&](const string& row) {
        const auto r = layout.Parse(row, default_person_id);
        pop->CreatePerson(r.person_id, r.age, r.profession, r.household, r.school, r.workplace,
                          r.community_weekend, r.community_weekday, r.household_cluster, r.collectivity);
        ++default_person_id;
    };
    if (!layout.HasHeader()) {
        add_person(line);
    }
    while (getline(popFile, line)) {
        add_person(line);
    }

    popFile.close();

    return pop;
}

std::string RemoveQuotes(const std::string& input) {
                if (!input.empty() && input.front() == '"' && input.back() == '"') {
                return input.substr(1, input.size() - 2);
                }
                return input;
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
                        if (typ != Id::OtherHouse && typ != Id::RestoCafe && typ != Id::OtherPlace && typ != Id::Transport) {
                        maxIds[typ] = max(maxIds[typ], p.GetPoolId(typ));
                        }

                }
        }
        // --------------------------------------------------------------
        // Initialize poolSys with empty ContactPools (even for Id=0).
        // --------------------------------------------------------------
        for (Id typ : IdList) {
                if (typ != Id::OtherHouse && typ != Id::RestoCafe && typ != Id::OtherPlace && typ != Id::Transport) {
                for (unsigned int i = 1; i < maxIds[typ] + 1; i++) {
                        pop->RefPoolSys().CreateContactPool(typ);
                }}
        }

        // --------------------------------------------------------------
        // Insert persons (pointers) in their contactpools. Having Id 0
        // means "not belonging pool of that type" (e.g. school/ workplace -
        // cannot belong to both, or e.g. out-of-work).
        //
        // Pools are uniquely identified by (type, subscript) and a Person
        // belongs, per type, to the pool with subscript p.GetPoolId(type).
        // Defensive measure: we have a pool for Id 0 and leave it empty.
        // --------------------------------------------------------------
        
        std::map<int, Person*> id_pointer_persons;
        
        for (auto& p : *pop) {
                unsigned int person_id = p.GetId();
                id_pointer_persons.emplace(person_id, &p);

                for (Id typ : IdList) {
                        if (typ != Id::OtherHouse && typ != Id::RestoCafe && typ != Id::OtherPlace && typ != Id::Transport) {
                        const auto poolId = p.GetPoolId(typ);
                        if (poolId > 0) {
                                pop->RefPoolSys().RefPools(typ)[poolId].AddMember(&p);
                        }
                }}
                
        }

        const auto allowed_subpools_communities = m_config.get<bool>("run.subpools_community_used",false);
        if (allowed_subpools_communities) {

        const auto fileName = m_config.get<string>("run.subpools_community_file");
        m_stride_logger->info("Building subpools from file {}.", fileName);
        const filesys::path filePath{fileName};
        if (!is_regular_file(filePath)) {
        throw runtime_error(string(__func__) + "> subpools community file " + filePath.string() + " not present.");
        }
    
                ifstream subpoolsCommunityFile;
                subpoolsCommunityFile.open(filePath.string());
                if (!subpoolsCommunityFile.is_open()) {
                throw runtime_error(string(__func__) + "> Error opening new community file " + filePath.string());
                }



                string line;
                getline(subpoolsCommunityFile, line); // step over file header
                auto headers   = Split(line, ";");

                int line_number = 1;

                while (getline(subpoolsCommunityFile, line)) {
                const auto values               = Split(line, ";");
                 
        const auto person_id = static_cast<unsigned int>(IntFromString(values[0]));
        const auto subpool_id = static_cast<unsigned int>(IntFromString(values[1]));
        const std::string& location = values[2];
        const auto day_week = static_cast<unsigned int>(IntFromString(values[3]));
        const auto duration = static_cast<unsigned int>(IntFromString(values[4]));
                          
                ContactType::Id typ = ToId(location);

                if (line_number < 5) {
                maxIds[typ] = subpool_id;
                        } 
                if (line_number == 5) {
                // --------------------------------------------------------------
        // Initialize poolSys with empty ContactPools (even for Id=0).
        // --------------------------------------------------------------
                for (Id typ : IdList) {
                        if (typ == Id::OtherHouse || typ == Id::RestoCafe || typ == Id::OtherPlace || typ == Id::Transport) {
                                for (unsigned int i = 1; i < maxIds[typ] + 1; i++) {
                                        
                                pop->RefPoolSys().CreateContactPool(typ);
                }}}

                }

                if (line_number > 4) {
                              
                Person* p=id_pointer_persons[person_id];
                               
                if (subpool_id > 0) {
                
                pop->RefPoolSys().RefPools(typ)[subpool_id].SetDayWeek(day_week);
                pop->RefPoolSys().RefPools(typ)[subpool_id].AddMember(p);
                }
                
                auto& record = pop->RefVenueAttendance().Ref(*p)[VenueAttendance::VenueIndex(typ)].at(day_week);
                record.pool_id  = subpool_id;
                record.duration = duration;

                }

                line_number++;
                
                }

                subpoolsCommunityFile.close();

                // Venue pools now have their final members and day: give them the
                // durations (contacts follow in ContactDivider).
                pop->RefVenueAttendance().CopyToPools(pop->RefPoolSys());
       
        }
        

       
m_stride_logger->info("Building population ready");
        


        return pop;
}

} // namespace stride
