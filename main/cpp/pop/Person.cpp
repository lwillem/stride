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
 * Header file for the Person class.
 */

#include "Person.h"

#include "contact/ContactType.h"
#include "pop/Age.h"

#include <stdexcept>

namespace stride {

using namespace std;
using namespace stride::ContactType;
using namespace stride::util;

void Person::UpdateEvents(unsigned int simDay)
{
    if (!m_events.empty() && m_events.top().GetTime() < simDay) {
        throw std::runtime_error("Person event scheduled in the past!");
    }

    // Update events
    while (!m_events.empty() && m_events.top().GetTime() == simDay) {
        const Event e = m_events.top();
        m_events.pop();

        //TODO: log the execution of events, 
        //make the logger available through a singleton? 
        EventType type = e.GetType();
        if (type == Person::EventType::StartIsolation) {
            m_isolated = true;
        } else if (type == Person::EventType::EndIsolation) {
            m_isolated = false;
        }
    }
}

void Person::ScheduleEvent(unsigned int simDay, const Event &event)
{
    m_events.push(event);
    if (simDay == event.GetTime())
        UpdateEvents(simDay);
} 

void Person::Isolate(unsigned int simDay, unsigned int from, unsigned int to)
{
    auto start = Person::Event(from, Person::EventType::StartIsolation);
    ScheduleEvent(simDay, start);
    auto end = Person::Event(to, Person::EventType::EndIsolation);
    ScheduleEvent(simDay, end);
}

//TODO: boolean args can be obtained from the calendar
void Person::Update(bool isRegularWeekday, bool isK12SchoolOff, bool isCollegeOff,
		bool isHouseholdClusteringAllowed,
        bool isIsolatedFromHousehold,
		util::RnHandler& rnHandler,
        const std::shared_ptr<Calendar> calendar)

{
        const unsigned int simDay = calendar->GetSimulationDay();

        UpdateEvents(simDay);

        // Update health and disease status
        m_health.Update();

        // by default: a person is at home (or in their collectivity)
        m_in_pools[Id::Household]          = true;
        m_in_pools[Id::Collectivity]       = true;


        // is household clustering allowed?
        m_in_pools[Id::HouseholdCluster]   = isHouseholdClusteringAllowed ? true : false;

        // Update presence in contact pools by type of day
        if (isRegularWeekday) {
        	m_in_pools[Id::Workplace]          = true;
			m_in_pools[Id::PrimaryCommunity]   = false;
			m_in_pools[Id::SecondaryCommunity] = true;
        } else{
            m_in_pools[Id::Workplace]          = false;
            m_in_pools[Id::PrimaryCommunity]   = true;
            m_in_pools[Id::SecondaryCommunity] = false;
        }

        // Update presence at school and college
        m_in_pools[Id::K12School] = isK12SchoolOff ? false : true;
        m_in_pools[Id::College]   = isCollegeOff   ? false : true;


        // Update presence in contact pools by health state
        if (m_health.IsSymptomatic()) {

        	// probability of staying home from school/work given symptoms
        	if(rnHandler.Binomial(m_health.GetSymptomaticCntReductionWorkSchool())){
        		m_in_pools[Id::K12School]          = false;
				m_in_pools[Id::College]            = false;
				m_in_pools[Id::Workplace]          = false;
        	}

            // probability of staying home from community pools given symptoms
        	if(rnHandler.Binomial(m_health.GetSymptomaticCntReductionCommunity())){
				m_in_pools[Id::PrimaryCommunity]   = false;
				m_in_pools[Id::SecondaryCommunity] = false;
        	}

        	// stay home from household cluster when symptomatic
			m_in_pools[Id::HouseholdCluster]   = false;

        }

        // Update presence in contact pools if person is in quarantine
        if(InIsolation()){
        	m_in_pools[Id::Household]          = !isIsolatedFromHousehold;  
        	m_in_pools[Id::K12School]          = false;
			m_in_pools[Id::College]            = false;
			m_in_pools[Id::Workplace]          = false;
        	m_in_pools[Id::PrimaryCommunity]   = false;
        	m_in_pools[Id::SecondaryCommunity] = false;
        	m_in_pools[Id::HouseholdCluster]   = false;
        	m_in_pools[Id::Collectivity]       = false;  //TODO: correct assumption?!
        }

} // Person::Update()

bool operator>(const Person::Event& lhs, const Person::Event& rhs)
{
    return lhs.GetTime() > rhs.GetTime();
}

} // namespace stride
