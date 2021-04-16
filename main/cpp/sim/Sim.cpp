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
 *  Copyright 2020, Kuylen E, Willem L, Broeckhove J
 */

/**
 * @file
 * Implementation for the Simulator class.
 */

#include "Sim.h"

#include "calendar/Calendar.h"
#include "contact/ContactType.h"
#include "contact/InfectorExec.h"
#include "disease/DiseaseSeeder.h"
#include "pop/Population.h"
#include "sim/SimBuilder.h"
#include "util/RunConfigManager.h"
#include "util/StringUtils.h"

#include <omp.h>
#include <utility>

namespace stride {

using namespace std;
using namespace stride::util;
using namespace EventLogMode;

Sim::Sim()
    : m_config(), m_event_log_mode(Id::None), m_num_threads(1U), m_track_index_case(false),
      m_calendar(nullptr), m_contact_profiles(), m_rn_handlers(), m_infector_default(),m_infector_tracing(),
      m_population(nullptr), m_rn_man(), m_transmission_profile(),
	  m_cnt_intensity_householdCluster(0),
      m_is_isolated_from_household(false),
	  m_public_health_agency(),m_universal_testing(),m_num_daily_imported_cases(0),m_hospitalisation_config()
{
}

std::shared_ptr<Sim> Sim::Create(const boost::property_tree::ptree& config, shared_ptr<Population> pop,
                                 util::RnMan rnMan)
{
        struct make_shared_enabler : public Sim
        {
                explicit make_shared_enabler() : Sim() {}
        };
        shared_ptr<Sim> sim = make_shared<make_shared_enabler>();
        SimBuilder(config).Build(sim, std::move(pop), std::move(rnMan));
        return sim;
}

void Sim::TimeStep()
{

        // Logic where you compute (on the basis of input/config for initial day or on the basis of
        // number of sick persons, duration of epidemic etc) what kind of DaysOff scheme you apply.
        const bool isRegularWeekday     = m_calendar->IsRegularWeekday();
        const bool isHouseholdClusteringAllowed    = m_calendar->IsHouseholdClusteringAllowed();

		// To be used in update of population & contact pools.
        Population& population    = *m_population;
        auto&       logger        = population.RefEventLogger();
        auto&       poolSys       = population.RefPoolSys();
        auto        eventLogger   = population.RefEventLogger();
        const auto  simDay        = m_calendar->GetSimulationDay();

        // Select infector, based on tracing
        const auto& infector      = m_public_health_agency.IsContactTracingActive(m_calendar) ? *m_infector_tracing : *m_infector_default;

        // set HouseholdCluster intensity
        double cnt_intensity_householdCluster = 0.0;
		if (isHouseholdClusteringAllowed && poolSys.RefPools(ContactType::Id::HouseholdCluster).size() > 1){
			cnt_intensity_householdCluster = m_cnt_intensity_householdCluster;
		}

        // Import infected cases into the population
        if(m_calendar->GetNumberOfImportedCases() > 0){
        	DiseaseSeeder(m_config, m_rn_man).ImportInfectedCases(m_population, m_calendar->GetNumberOfImportedCases(), simDay, m_transmission_profile, m_rn_handlers[0]);
            logger->info("[IMPORT-CASES] sim_day={} count={}", simDay, m_calendar->GetNumberOfImportedCases());        	
        }

#pragma omp parallel num_threads(m_num_threads)
        {
        	const auto thread_num = static_cast<unsigned int>(omp_get_thread_num());
			// Update health status and presence/absence in contact pools
			// depending on health status, work/school day and whether
			// we want to track index cases without adaptive behavior
#pragma omp for schedule(static)
			for (size_t i = 0; i < population.size(); ++i) {

				// adjust K12SchoolOff boolean to school type for individual 'i'
//				bool isK12SchoolOff = m_public_health_agency.IsK12SchoolOff(population[i].GetAge(),
//						m_calendar->IsSchoolClosed(0), //isPreSchoolOff,
//						m_calendar->IsSchoolClosed(6), //isPrimarySchoolOff,
//						m_calendar->IsSchoolClosed(11), //isSecondarySchoolOff,
//						m_calendar->IsSchoolClosed(20)); //isCollegeOff);

				unsigned int school_id = population[i].GetPoolId(ContactType::Id::K12School);
				unsigned int school_age = population[i].GetAge();
				if(school_id>0){
					school_age = poolSys.RefPools(ContactType::Id::K12School)[school_id].GetMinAge();
				}
				bool isK12SchoolOff = m_calendar->IsSchoolClosed(school_age);
				bool isCollegeOff   = m_calendar->IsSchoolClosed(population[i].GetAge());
				// update health and presence at different contact pools
				population[i].Update(isRegularWeekday, isK12SchoolOff, isCollegeOff,
						isHouseholdClusteringAllowed,
						m_is_isolated_from_household,
                        m_rn_handlers[thread_num], 
                        m_calendar);
			}
        }// end pragma openMP

		 // Perform contact tracing (if activated)
		 m_public_health_agency.PerformContactTracing(m_population, m_rn_handlers[0], m_calendar);

		 // Perform universal testing 
	     m_universal_testing.PerformUniversalTesting(m_population, m_rn_handlers[0], m_calendar,m_public_health_agency);

#pragma omp parallel num_threads(m_num_threads)
        {
		    const auto thread_num = static_cast<unsigned int>(omp_get_thread_num());
			// Infector updates individuals for contacts & transmission within each pool.
			// Skip pools with id = 0, because it means Not Applicable.
			for (auto typ : ContactType::IdList) {
					if ((typ == ContactType::Id::Workplace && !isRegularWeekday) ||
						(typ == ContactType::Id::K12School && !isRegularWeekday) ||
						(typ == ContactType::Id::College && !isRegularWeekday) ||
						(typ == ContactType::Id::HouseholdCluster && !isHouseholdClusteringAllowed)) {
							continue;
					}
#pragma omp for schedule(static)
					for (size_t i = 1; i < poolSys.RefPools(typ).size(); i++) { // NOLINT
							infector(poolSys.RefPools(typ)[i], m_contact_profiles[typ], m_transmission_profile,
									 m_rn_handlers[thread_num], simDay, eventLogger,
									 m_population,cnt_intensity_householdCluster,m_calendar);
					}
			}
        } // end pragma openMP

        m_population->RefEventLogger()->flush();
        m_calendar->AdvanceDay();
}

} // namespace stride
