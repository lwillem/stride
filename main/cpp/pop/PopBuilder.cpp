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
        auto headers   = Split(line, ";");

        while (getline(popFile, line)) {
                const auto values               = Split(line, ";");
                const auto age                  = FromString<unsigned int>(values[0]);
                const auto person_id            = FromString<unsigned int>(values[1]);
                const auto householdId          = FromString<unsigned int>(values[2]);
                auto schoolId                   = FromString<unsigned int>(values[3]);
                const auto workId               = FromString<unsigned int>(values[4]);
                const auto primaryCommunityId   = FromString<unsigned int>(values[5]);
                const auto secondaryCommunityId = FromString<unsigned int>(values[6]);

                unsigned int householdClusterId = 0;
                if(values.size() == 8 && Trim(ToString(headers[7]),ToString('"')) == "household_cluster_id"){
                	householdClusterId = FromString<unsigned int>(values[7]);
                }

                unsigned int collectivityId = 0;

				if(values.size() == 8 && Trim(ToString(headers[7]),ToString('"')) == "collectivity_id"){
					collectivityId = FromString<unsigned int>(values[7]);
				}

                //TODO: rename school types to current approach
                unsigned int collegeId = 0;
                if(schoolId != 0 && age >= age_break_school_types && age < 23){
                	collegeId = schoolId;
                	schoolId = 0;
                }

                pop->CreatePerson(person_id, age, householdId, schoolId, collegeId, workId, primaryCommunityId,
                                  secondaryCommunityId, householdClusterId, collectivityId);
                ;
        }

        popFile.close();

        return pop;
}

shared_ptr<Population> PopBuilder::MakePersonsOpt(shared_ptr<Population> pop)
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
    auto headers   = Split(line, ";");
    //
    bool has_extra_column = headers.size() == 8;
    string extra_id = "";
    if (has_extra_column) { extra_id = Trim(ToString(headers[7]),ToString('"')); }
    bool household_cluster_id = extra_id == "household_cluster_id";
    bool collectivity_id = extra_id == "collectivity_id";
    const unsigned int defaultHouseholdClusterId = 0;
    const unsigned int defaultCollectivityId = 0;

    // Read lines from file
    

    while (getline(popFile, line)) {
        const auto values               = Split(line, ";");
        const auto age                  = static_cast<unsigned int>(IntFromString(values[0]));
        const auto person_id            = static_cast<unsigned int>(IntFromString(values[1]));
        const auto householdId          = static_cast<unsigned int>(IntFromString(values[2]));
        auto schoolId                   = static_cast<unsigned int>(IntFromString(values[3]));
        const auto workId               = static_cast<unsigned int>(IntFromString(values[4]));
        const auto primaryCommunityId   = static_cast<unsigned int>(IntFromString(values[5]));
        const auto secondaryCommunityId = static_cast<unsigned int>(IntFromString(values[6]));
       
        unsigned int householdClusterId = defaultHouseholdClusterId;
        unsigned int collectivityId = defaultCollectivityId;
        if (values.size() == 8) {
            if (household_cluster_id) {
                householdClusterId = static_cast<unsigned int>(IntFromString(values[7]));
            } else if (collectivity_id) {
                collectivityId = static_cast<unsigned int>(IntFromString(values[7]));
            }
        }

        //TODO: rename school types to current approach
        unsigned int collegeId = 0;
        if(schoolId != 0 && age >= age_break_school_types && age < 23){
            collegeId = schoolId;
            schoolId = 0;
        }

        if(person_id == 1 ) {
                  std::cout << "line " << line << std::endl;
                  
                }

        pop->CreatePerson(person_id, age, householdId, schoolId, collegeId, workId, primaryCommunityId,
                          secondaryCommunityId, householdClusterId, collectivityId);
        ;
    }


    popFile.close();

    std::cout << "popFile closed" << std::endl;

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
//        MakePersons(pop);  // Toggle between these two if preferring old version
        MakePersonsOpt(pop);

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
                        pop->RefPoolSys().CreateContactPool(typ, 7);
                }}
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


        const auto allowed_subpools_communities = m_config.get<bool>("run.subpools_community_used");
        const auto fileName = m_config.get<string>("run.subpools_community_file");
        m_stride_logger->info("Building subpools from file {}.", fileName);
        const auto use_install_dirs = m_config.get<bool>("run.use_install_dirs");
        const auto filePath         = (use_install_dirs) ? FileSys::GetDataDir() /= fileName : filesys::path(fileName);
        if (!is_regular_file(filePath)) {
        throw runtime_error(string(__func__) + "> subpools community file " + filePath.string() + " not present.");
        }

        if (allowed_subpools_communities) {
    
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

        // printen tijdens inlezen bestand
                  if (subpool_id > 0)   {          
                cout << "Processing line " << line_number << ": " << line << endl; // Print lijnnummer en inhoud van de lijn
                  }
           
                ContactType::Id typ = ToId(location);

                // unsigned int last_typ_pool = pop->RefPoolSys().currentPoolIds(typ);
   

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
                                pop->RefPoolSys().CreateContactPool(typ, day_week);
                }}}

                }

                if (line_number > 4) {
              
                //if (subpool_id > last_typ_pool) {
                //       pop->RefPoolSys().CreateContactPool(typ,day_week);
                //unsigned int new_last_typ_pool = pop->RefPoolSys().currentPoolIds(typ);

                //}
           
                Person* p=id_pointer_persons[person_id];
                               
                if (subpool_id > 0) {
                pop->RefPoolSys().RefPools(typ)[subpool_id].AddMember(p);
                }
                
                Person person = *p;
                person.PoolIds(typ)[day_week] = subpool_id;
                person.PoolDurations(typ)[day_week] = duration;

                }

                line_number++;
                
                }

                subpoolsCommunityFile.close();
       
        }
        

       
m_stride_logger->info("Building population ready");
        


        return pop;
}

} // namespace stride
