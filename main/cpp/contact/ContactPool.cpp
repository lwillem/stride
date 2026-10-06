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
 *  Copyright 2019, Kuylen E, Willem L, Broeckhove J
 */

/**
 * @file
 * Implementation for the core ContcatPool class.
 */

#include "ContactPool.h"

#include "pop/Age.h"
#include "pop/Person.h"

#include <algorithm>

namespace stride {

using namespace std;

ContactPool::ContactPool(unsigned int poolId, ContactType::Id type)
    : m_index_immune(0), m_pool_id(poolId), m_pool_type(type), m_members(), m_min_age(110), m_ventilation(0), m_venue_non_complier(false), m_day_week(7), m_air_mass(1.0),
      m_duration(0), m_member_durations(), m_member_contacts()
{
}

void ContactPool::AddMember(Person* p)
{
        m_members.emplace_back(p);
        m_index_immune++;

        if(p->GetAge()<m_min_age){
        	m_min_age = p->GetAge();
        }
}

void ContactPool::SetMemberAttendance(size_t i, unsigned int duration, unsigned int contacts)
{
        if (m_member_durations.size() != m_members.size()) {
                m_member_durations.assign(m_members.size(), 0U);
                m_member_contacts.assign(m_members.size(), 0U);
        }
        m_member_durations[i] = duration;
        m_member_contacts[i]  = contacts;
}

void ContactPool::SwapMembers(size_t i, size_t j)
{
        swap(m_members[i], m_members[j]);
        if (!m_member_durations.empty()) {
                swap(m_member_durations[i], m_member_durations[j]);
                swap(m_member_contacts[i], m_member_contacts[j]);
        }
}

unsigned int ContactPool::GetInfectedCount() const
{
        unsigned int infected = 0;

        for (stride::Person* person : m_members) {
                if (person->GetHealth().IsInfected()) {
                        infected++;
                }
        }
        return infected;
}

std::tuple<bool, unsigned int> ContactPool::SortMembers()
{
        bool         infectious_cases = false;
        unsigned int num_cases        = 0;

        for (size_t i_member = 0; i_member < m_index_immune; i_member++) {
                // if immune, move to back
                if (m_members[i_member]->IsImmune()) {
                        bool         swapped   = false;
                        unsigned int new_place = m_index_immune - 1;
                        m_index_immune--;
                        while (!swapped && new_place > i_member) {
                                if (m_members[new_place]->IsImmune()) {
                                        m_index_immune--;
                                        new_place--;
                                } else {
                                        SwapMembers(i_member, new_place);
                                        swapped = true;
                                }
                        }
                }
                // else, if not susceptible, move to front
                else if (!m_members[i_member]->GetHealth().IsSusceptible()) {
                        if (!infectious_cases && m_members[i_member]->GetHealth().IsInfectious()) {
                                infectious_cases = true;
                        }
                        if (i_member > num_cases) {
                                SwapMembers(i_member, num_cases);
                        }
                        num_cases++;
                }
        }
        return std::make_tuple(infectious_cases, num_cases);
}

} // namespace stride
