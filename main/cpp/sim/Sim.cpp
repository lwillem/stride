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
	  m_run_simplified(false),
      m_calendar(nullptr), m_contact_profiles(), m_rn_handlers(), m_infector_default(),m_infector_tracing(),
      m_population(nullptr), m_rn_man(), m_transmission_profile(),
      m_is_isolated_from_household(false),
	  m_public_health_agency(),m_hospitalisation_config()
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
        // Define the type of day
        const bool isRegularWeekday                = m_calendar->IsRegularWeekday();

		// To be used in population & contact pool update
        Population& population    = *m_population;
        auto&       logger        = population.RefEventLogger();
        auto&       poolSys       = population.RefPoolSys();
        auto        eventLogger   = population.RefEventLogger();
        const auto  simDay        = m_calendar->GetSimulationDay();

        // Select infector, based on tracing
        const auto& infector      = m_public_health_agency.IsContactTracingActive(m_calendar) ? *m_infector_tracing : *m_infector_default;

        // Get household clustering intensity
        double cnt_intensity_householdCluster = m_calendar->GetHouseholdClusteringLevel();

        // Update Health before introducing new cases (infected on simDay)
#pragma omp parallel num_threads(m_num_threads)
        {
#pragma omp for schedule(static)
        	for (size_t i = 0; i < population.size(); ++i) {
        		population[i].UpdateHealth();

			}
		} // end pragma openMP

        // Import infected cases into the population
        if(m_calendar->GetNumberOfImportedCases() > 0){
        	DiseaseSeeder(m_config, m_rn_man).ImportInfectedCases(m_population, m_calendar->GetNumberOfImportedCases(), simDay, m_transmission_profile, m_rn_handlers[0]);
            logger->info("[IMPORT-CASES] sim_day={} count={}", simDay, m_calendar->GetNumberOfImportedCases());        	
        }

#pragma omp parallel num_threads(m_num_threads)
        {
        	const auto thread_num = static_cast<unsigned int>(omp_get_thread_num());
			// Update presence/absence in contact pools depending on health status
#pragma omp for schedule(static)
			for (size_t i = 0; i < population.size(); ++i) {

				// update health-related presence at different contact pools
				population[i].UpdatePresence(m_is_isolated_from_household,
                        m_rn_handlers[thread_num], 
                        simDay, m_run_simplified);
			}
        }// end pragma openMP

		 // Perform contact tracing (if activated)
		 m_public_health_agency.PerformContactTracing(m_population, m_rn_handlers, m_calendar);

		 // Process social contact behaviour and transmission dynamics
#pragma omp parallel num_threads(m_num_threads)
        {
		    const auto thread_num = static_cast<unsigned int>(omp_get_thread_num());
			// Infector updates individuals for contacts & transmission within each pool.
		    // Skip Communities, Workplaces, Schools or HouseholdClusters is possible.
			for (auto typ : ContactType::IdList) {
					if ((typ == ContactType::Id::PrimaryCommunity && isRegularWeekday) ||
						(typ == ContactType::Id::SecondaryCommunity && !isRegularWeekday) ||
						(typ == ContactType::Id::Workplace && !isRegularWeekday) ||
						(typ == ContactType::Id::School && !isRegularWeekday) ||
						(typ == ContactType::Id::HouseholdCluster && cnt_intensity_householdCluster==0)) {
							continue;
					}
#pragma omp for schedule(static)
					// Skip pools with id = 0, because it means Not Applicable.
					for (size_t i = 1; i < poolSys.RefPools(typ).size(); i++) { // NOLINT

						// enable ContactPool specific physical distancing
						double typ_distancing_factor = m_calendar->GetDistancingFactor(poolSys.RefPools(typ)[i]);

						infector(poolSys.RefPools(typ)[i], m_contact_profiles[typ], m_transmission_profile,
								 m_rn_handlers[thread_num], simDay, eventLogger,
								 m_population, cnt_intensity_householdCluster, typ_distancing_factor);
					}
			}
        } // end pragma openMP

//        // log prevalence? (time consuming!)
//        m_population->LogPrevalence(simDay);

        // flush event logger and advance one day
        m_population->RefEventLogger()->flush();
        m_calendar->AdvanceDay();
}

} // namespace stride
