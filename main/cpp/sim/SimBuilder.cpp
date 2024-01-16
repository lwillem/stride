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
#include "contact/VentilationHeterogeneitySeeder.h"
#include "disease/DiseaseSeeder.h"
#include "disease/HealthSeeder.h"
#include "disease/ImmunitySeeder.h"
#include "disease/PublicHealthAgency.h"
#include "pop/SurveySeeder.h"
#include "sim/Sim.h"
#include "util/FileSys.h"
#include "util/RnMan.h"

namespace stride {

using namespace boost::property_tree;
using namespace std;
using namespace stride::util;
using namespace ContactType;

SimBuilder::SimBuilder(const ptree& config) : m_config(config) {}

shared_ptr<Sim> SimBuilder::Build(shared_ptr<Sim> sim, shared_ptr<Population> pop, RnMan rnMan)
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
        sim->m_airborne_transmission          = m_config.get<bool>("run.airborne_transmission", false);
        // TO DO!!! test only airborne transmission if subpools community!!!
        sim->m_num_threads                   = m_config.get<unsigned int>("run.num_threads");
        unsigned int num_days                = m_config.get<unsigned short>("run.num_days");
        sim->m_calendar                      = make_shared<Calendar>(m_config,num_days);
        sim->m_event_log_mode                = EventLogMode::ToMode(m_config.get<string>("run.event_log_level", "None"));
        sim->m_rn_man                        = std::move(rnMan);

        std::cout << "Contact handlers" << std::endl;
        // --------------------------------------------------------------
        // Contact handlers, each with generator bound to different
        // random engine stream) and infector.
        // --------------------------------------------------------------
        for (unsigned int i = 0; i < sim->m_num_threads; i++) {
                auto gen = sim->m_rn_man.GetUniform01Generator(i);
                sim->m_rn_handlers.emplace_back(util::RnHandler(gen));
        }
        const auto& select = make_tuple(sim->m_event_log_mode, sim->m_track_index_case);
        sim->m_infector_default    = InfectorMap().at(select);

        // additional infector if logmode is Tracing
        if(m_config.get<string>("run.event_log_level", "None") == "ContactTracing"){
        	const auto& select_tracing  = make_tuple(EventLogMode::Id::All, sim->m_track_index_case);
        	sim->m_infector_tracing    = InfectorMap().at(select_tracing);
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
        
         std::cout << "Initialize the airborne transmission profile." << std::endl;
        // --------------------------------------------------------------
        // Initialize the airborne transmission profile
        // --------------------------------------------------------------
        const auto diseasePt = ReadDiseaseAirborneTransmissionPtree();
        sim->m_airborne_transmission_profile.Initialize(m_config, diseasePt);

        std::cout << "Seed the population with health data." << std::endl;
        // --------------------------------------------------------------
        // Seed the population with health data.
        // --------------------------------------------------------------
        HealthSeeder(diseasePt).Seed(sim->m_population, sim->m_transmission_profile, sim->m_rn_handlers);

        std::cout << "Seed population with immunity: naturally or vaccine-induced." << std::endl;
        // --------------------------------------------------------------
	// Seed population with immunity: naturally or vaccine-induced.
	// --------------------------------------------------------------
        ImmunitySeeder(m_config, sim->m_rn_man).Seed(sim->m_population);

        std::cout << "Seed population with infection." << std::endl;
        // --------------------------------------------------------------
        // Seed population with infection.
        // --------------------------------------------------------------
        DiseaseSeeder(m_config, sim->m_rn_man).Seed(sim->m_population, sim->m_transmission_profile, sim->m_rn_handlers[0]);
        sim->m_num_daily_imported_cases = m_config.get<double>("run.num_daily_imported_cases",0);

        std::cout << "Set Universal Testing " << std::endl;
        // --------------------------------------------------------------
	// Set Universal Testing 
        // --------------------------------------------------------------
        sim->m_universal_testing.Initialize(m_config);
        
        std::cout << "Set Public Health Agency" << std::endl;
        // --------------------------------------------------------------
        // Set Public Health Agency
        // --------------------------------------------------------------
        sim->m_public_health_agency.Initialize(m_config);
		sim->m_cnt_intensity_householdCluster       = m_config.get<double>("run.cnt_intensity_householdCluster",0);
		sim->m_is_isolated_from_household           = m_config.get<bool>("run.is_isolated_from_household",false);

        std::cout << "Seed population with survey participants." << std::endl;
        // --------------------------------------------------------------
        // Seed population with survey participants.
        // --------------------------------------------------------------
        SurveySeeder(m_config, sim->m_rn_man).Seed(sim->m_population);

        std::cout << "Seed population with non-compliant individuals." << std::endl;
        // --------------------------------------------------------------
        // Seed population with non-compliant individuals.
        // --------------------------------------------------------------
        ContactHeterogeneitySeeder(m_config, sim->m_rn_man).Seed(sim->m_population);

        std::cout << "Calculate contacts based on age contact profile and duration in location" << std::endl;
        //---------------------------------------------------------------
        // Calculate contacts based on age contact profile and duration in location
        //---------------------------------------------------------------
        if (sim->m_subpools_community){
                ContactDivider(m_config, sim->m_rn_man).Divide(sim->m_population, sim->m_contact_profiles);
        };

        std::cout << "Fill in the ventilation in the contactPoolSys." << std::endl;
        // --------------------------------------------------------------
        // Fill in the ventilation in the contactPoolSys.
        // --------------------------------------------------------------
        VentilationHeterogeneitySeeder(m_config, sim->m_rn_man).Ventilation(sim->m_population);

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

ptree SimBuilder::ReadDiseaseAirborneTransmissionPtree()
{
        const auto fn = m_config.get<string>("run.disease_airborne_transmission_file");
        const auto fp = m_config.get<bool>("run.use_install_dirs") ? FileSys::GetDataDir() /= fn : filesys::path(fn);
        return FileSys::ReadPtreeFile(fp);
}

} // namespace stride
