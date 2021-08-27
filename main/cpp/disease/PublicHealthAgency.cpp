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
 *  Copyright 2020, Willem L
 */

/**
 * @file
 * Implementation for the PublicHealthAgency class.
 */

#include "PublicHealthAgency.h"

#include "calendar/Calendar.h"
#include "pop/Population.h"
#include "util/FileSys.h"
#include "util/LogUtils.h"
#include "util/StringUtils.h"

#include <boost/property_tree/ptree.hpp>
#include <algorithm>

namespace stride {

using namespace boost::property_tree;
using namespace stride::ContactType;
using namespace stride::util;
using namespace std;

// Default constructor
PublicHealthAgency::PublicHealthAgency(): m_detection_probability(0),
		m_tracing_efficiency_household(0),m_tracing_efficiency_other(0),m_case_finding_capacity(0),m_delay_isolation_index(0),m_delay_contact_tracing(0),
		m_test_false_negative(0)
	{}

void PublicHealthAgency::Initialize(const ptree& config){
	m_detection_probability       = config.get<double>("run.detection_probability",0);

	m_tracing_efficiency_household = config.get<double>("run.tracing_efficiency_household",0);
	m_tracing_efficiency_other     = config.get<double>("run.tracing_efficiency_other",0);

	m_case_finding_capacity  = config.get<unsigned int>("run.case_finding_capacity",0);

	m_delay_isolation_index          = config.get<unsigned int>("run.delay_isolation_index",3);
	m_delay_contact_tracing  = config.get<unsigned int>("run.delay_contact_tracing",1);
	m_test_false_negative    = config.get<double>("run.test_false_negative",0.3);

	// account for false negative tests
	m_detection_probability        *= (1.0 - m_test_false_negative);
	m_tracing_efficiency_household  *= (1.0 - m_test_false_negative);
	m_tracing_efficiency_other      *= (1.0 - m_test_false_negative);
}

bool PublicHealthAgency::IsContactTracingActive(const std::shared_ptr<Calendar> calendar) const {
	return (m_detection_probability > 0) && calendar->IsContactTracingActivated();
}


void PublicHealthAgency::PerformContactTracing(std::shared_ptr<Population> pop, RnHandler& rnHandler,
												const std::shared_ptr<Calendar> calendar)
{

	// if contact tracing not active, stop
	if (!IsContactTracingActive(calendar)) {
			return;
	}

	//cout << m_detection_probability << " -- "<< m_tracing_efficiency_household << " -- "<< m_tracing_efficiency_other << " ** " << m_case_finding_capacity << endl;


	/// Mark index cases for track&trace
	for (auto& p_case : *pop) {

		if(p_case.GetHealth().NumberDaysInfected(1) &&
				rnHandler.Binomial(m_detection_probability)) {
			p_case.SetTracingIndexCase();
		}

	}

	/// Set counter for index cases
	unsigned int num_index_cases = 0;

	/// Loop over the population to find index cases on day X after symptom onset
	for (auto& p_case : *pop) {

		if (p_case.IsTracingIndexCase() && p_case.GetHealth().NumberDaysSymptomatic(m_delay_isolation_index)	) {
        Trace(p_case, pop, rnHandler, calendar);

        // update index case counter, and terminate if quota is reached
        num_index_cases++;
        if(num_index_cases >= m_case_finding_capacity){
            return;
        }
       } 
	}
}

//TODO: rename IsolateAndTrace()
void PublicHealthAgency::Trace(Person& p_case, 
        std::shared_ptr<Population> pop, 
		RnHandler& rnHandler,
        const std::shared_ptr<Calendar> calendar)
{
	auto& logger       = pop->RefEventLogger();
	const auto  simDay = calendar->GetSimulationDay();

			// Set index case in quarantine.
		    // As this individual tested positive, he/she is isolated for 7 days.
		     	unsigned int start = simDay + 1; //start tomorrow
			p_case.Isolate(simDay, start, start+7);

			// counter for number of contacts tested
			unsigned int num_contacts_tested = 0;

			// get contact register and iterator
			vector<Person*> cnt_register = p_case.GetContactRegister();
			vector<Person*>::iterator ip;

			// cout << p_case.GetId() << ": size (orig) " <<  cnt_register.size() << endl;

			// Sorting the vector, get the unique elements, and remove the others
			std::sort(cnt_register.begin(), cnt_register.end());
			ip = std::unique(cnt_register.begin(), cnt_register.end());
			cnt_register.resize(std::distance(cnt_register.begin(), ip));

			// cout << p_case.GetId() << ": size (unique) " <<  cnt_register.size() << endl;

			// loop over contact register
			for (ip = cnt_register.begin(); ip != cnt_register.end(); ++ip) {


				Person* p_contact = *ip;

				// combine contact tracing efficiency and false negative rate
				double tracing_efficiency = m_tracing_efficiency_other;

				// set default poolType as "other
				std::string poolTypeString = "School";

				// if contact is part of same household, change tracing efficiency and poolType
				if(p_contact->GetPoolId(Id::Household) == p_case.GetPoolId(Id::Household)){
					poolTypeString    = ToString(Id::Household);
					tracing_efficiency = m_tracing_efficiency_household;
				}

				// if contact is part of same community, change poolType
				if(p_contact->GetPoolId(Id::PrimaryCommunity) == p_case.GetPoolId(Id::PrimaryCommunity)||
						p_contact->GetPoolId(Id::SecondaryCommunity) == p_case.GetPoolId(Id::SecondaryCommunity)){
					poolTypeString    = "Community";
				}

				// if contact is part of same workplace, change poolType
				if(p_case.GetPoolId(Id::Workplace) != 0 &&
						p_contact->GetPoolId(Id::Workplace) == p_case.GetPoolId(Id::Workplace)){
					poolTypeString    = "Workplace";
				}

				if(rnHandler.Binomial(tracing_efficiency)){

					if(p_contact->GetHealth().IsInfected()){
						// start isolation over X days
					    unsigned int start = simDay + m_delay_contact_tracing;
						p_contact->Isolate(simDay, start, start + 7);

						// add to log (TODO: check log_level)
						logger->info("[TRACE] {} {} {} {} {} {} {} {} {} {}",
								p_contact->GetId(), p_contact->GetAge(),
								p_contact->GetHealth().IsInfected(),
								p_contact->GetHealth().IsSymptomatic(),
								poolTypeString, p_case.GetId(),
								p_case.GetAge(), simDay, -1, -1);
					}
					// increment contact counter
					num_contacts_tested++;
				}
			}

			// Log index case
			// TODO: check log_level
			logger->info("[TRACE] {} {} {} {} {} {} {} {} {} {}",
						 p_case.GetId(), p_case.GetAge(),
						 p_case.GetHealth().IsInfected(),
						 p_case.GetHealth().IsSymptomatic(),
						 "Index", p_case.GetId(),
						 p_case.GetAge(), simDay, cnt_register.size(), num_contacts_tested);
}

} // namespace stride
