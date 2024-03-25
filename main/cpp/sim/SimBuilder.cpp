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
 */

/**
 * @file
 * Implementation for the SimBuilder class.
 */

#include "SimBuilder.h"

#include "contact/ContactType.h"
#include "contact/InfectorMap.h"
#include "contact/ContactDivider.h"
#include "contact/ContactHeterogeneitySeeder.h"
#include "contact/PoolCharacteristicsSeeder.h"
#include "health/DiseaseSeeder.h"
#include "health/HealthSeeder.h"
#include "health/ImmunitySeeder.h"
#include "healthcare/PublicHealthAgency.h"
#include "pop/SurveyManager.h"
#include "sim/Sim.h"
#include "util/StringUtils.h"
#include "util/FileSys.h"
#include "util/RnMan.h"

namespace stride {

using namespace boost::property_tree;
using namespace std;
using namespace stride::util;
using namespace ContactType;

SimBuilder::SimBuilder(const ptree& config) : m_config(config) {}

shared_ptr<Sim> SimBuilder::Build(shared_ptr<Sim> sim, shared_ptr<Population> pop)
{
        std::cout << "Read config info and setup random number manager" << std::endl;
        // --------------------------------------------------------------
        // Read config info and setup random number manager
        // --------------------------------------------------------------
        sim->m_config                        = m_config;
        sim->m_population                    = std::move(pop);
        sim->m_track_index_case              = m_config.get<bool>("run.track_index_case");
        sim->m_run_simplified                = m_config.get<bool>("run.run_simplified", false);
        sim->m_subpools_community            = m_config.get<bool>("run.subpools_community_used", false);
        sim->m_airborne_transmission         = m_config.get<bool>("run.airborne_transmission", false);
        // TO DO!!! test only airborne transmission if subpools community!!!
        sim->m_num_threads                   = m_config.get<unsigned int>("run.num_threads");
        unsigned int num_days                = m_config.get<unsigned short>("run.num_days");
        sim->m_calendar                      = make_shared<Calendar>(m_config,num_days);
        sim->m_event_log_mode                = EventLogMode::ToMode(m_config.get<string>("run.event_log_level", "None"));
        sim->m_rn_man_ptr 				     = std::make_shared<util::RnMan>(m_config.get<unsigned long>("run.rng_seed", 0U),
													                    	m_config.get<unsigned int>("run.num_threads"));


		// --------------------------------------------------------------
        // Select infector template(s) based on configuration.
        // --------------------------------------------------------------
        const auto& select = make_tuple(sim->m_event_log_mode, sim->m_track_index_case);
        sim->m_infector_default    = InfectorMap().at(select);

        // additional infector if logmode is ContactTracing or Participants
        if(m_config.get<string>("run.event_log_level", "None") == "ContactTracing"){
        	const auto& select_tracing  = make_tuple(EventLogMode::ToMode("ContactTracing"), sim->m_track_index_case);
        	sim->m_infector_tracing    = InfectorMap().at(select_tracing);
        } else if(m_config.get<string>("run.event_log_level", "None") == "Participants"){
        	sim->m_infector_tracing    = sim->m_infector_default;
        	const auto& select_default = make_tuple(EventLogMode::ToMode("Transmissions"), sim->m_track_index_case);
        	sim->m_infector_default    = InfectorMap().at(select_default);
        } else{
        	sim->m_infector_tracing    = InfectorMap().at(select);
        }

        std::cout << "Initialize the age-related contact profiles." << std::endl;
        // --------------------------------------------------------------
        // Initialize the age-related contact profiles.
        // --------------------------------------------------------------
        const auto ageContactPt = ReadAgeContactPtree();
        for (Id typ : IdList) {
                if (typ != Id::OtherHouse && typ != Id::RestoCafe && typ != Id::OtherPlace && typ != Id::Transport){
                sim->m_contact_profiles[typ] = AgeContactProfile(typ, ageContactPt);
                }
        }

        std::cout << "Initialize the transmission profile (fixes rates)." << std::endl;
        // --------------------------------------------------------------
        // Initialize the transmission profile (fixes rates).
        // --------------------------------------------------------------
        const auto diseasePt = ReadDiseasePtree();
        sim->m_transmission_profile.Initialize(m_config, diseasePt);
        
        std::cout << "Seed the population with health data." << std::endl;
        // --------------------------------------------------------------
        // Seed the population with health data (incl. hospital admission)
        // --------------------------------------------------------------
        HealthSeeder(m_config, diseasePt).Seed(sim->m_population, sim->m_transmission_profile, sim->m_rn_man_ptr);

        std::cout << "Seed population with immunity: naturally or vaccine-induced." << std::endl;
        // --------------------------------------------------------------
		// Seed population with immunity: naturally or vaccine-induced.
		// --------------------------------------------------------------
        ImmunitySeeder(m_config, sim->m_rn_man_ptr).Seed(sim->m_population);

        // --------------------------------------------------------------
        // Register infected seeds.
        // --------------------------------------------------------------
	sim->GetCalendar()->RegisterInfectedSeeds(m_config.get<unsigned int>("run.num_infected_seeds",0));

        std::cout << "Set Universal Testing " << std::endl;
        // --------------------------------------------------------------
		// Set Public Health Agency
		// --------------------------------------------------------------
        sim->m_public_health_agency.Initialize(m_config);
	sim->m_is_isolated_from_household = m_config.get<bool>("run.is_isolated_from_household",false);

        std::cout << "Seed population with survey participants." << std::endl;
        // --------------------------------------------------------------
        // Seed population with survey participants.
        // --------------------------------------------------------------
	sim->m_survey_manager = make_shared<SurveyManager>(sim->m_population, m_config, sim->m_rn_man_ptr);

        std::cout << "Seed population with non-compliant individuals." << std::endl;
        // --------------------------------------------------------------
        // Seed heterogeniety in social contact behaviour.
        // --------------------------------------------------------------
        ContactHeterogeneitySeeder(m_config, sim->m_rn_man_ptr).Seed(sim->m_population);

        std::cout << "Calculate contacts based on age contact profile and duration in location" << std::endl;
        //---------------------------------------------------------------
        // Calculate contacts based on age contact profile and duration in location
        //---------------------------------------------------------------
        if (sim->m_subpools_community){
                ContactDivider(m_config, sim->m_rn_man_ptr).Divide(sim->m_population, sim->m_contact_profiles);
        };

        std::cout << "Fill in characteristics in the contactPoolSys for airborne transmission." << std::endl;
        // --------------------------------------------------------------
        // Fill in characteristics in the contactPoolSys for airborne transmission
        // --------------------------------------------------------------
        const auto poolCharacteristicsPt = ReadPoolCharacteristicsPtree();
        PoolCharacteristicsSeeder(m_config, sim->m_rn_man_ptr).Seed(sim->m_population, poolCharacteristicsPt);

        // --------------------------------------------------------------
        // Done.
        // --------------------------------------------------------------
        return sim;
}

ptree SimBuilder::ReadAgeContactPtree()
{
        const auto fn = m_config.get<string>("run.age_contact_matrix_file", "contact_matrix.xml");
        const auto fp = m_config.get<bool>("run.use_install_dirs") ? FileSys::GetDataDir() /= fn : filesys::path(fn);
        return FileSys::ReadPtreeFile(fp);
}

ptree SimBuilder::ReadDiseasePtree()
{
        const auto fn = m_config.get<string>("run.disease_config_file");
        const auto fp = m_config.get<bool>("run.use_install_dirs") ? FileSys::GetDataDir() /= fn : filesys::path(fn);
        return FileSys::ReadPtreeFile(fp);
}

ptree SimBuilder::ReadPoolCharacteristicsPtree()
{
        const auto fn = m_config.get<string>("run.pool_characteristics_file");
        const auto fp = m_config.get<bool>("run.use_install_dirs") ? FileSys::GetDataDir() /= fn : filesys::path(fn);
        return FileSys::ReadPtreeFile(fp);
}


} // namespace stride
