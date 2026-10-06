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
 *  Copyright 2026, Kuylen E, Willem L, Broeckhove J
 */

/**
 * @file
 * Header file for the VenueAttendance class.
 */

#pragma once

#include "contact/ContactPoolSys.h"
#include "contact/ContactType.h"
#include "pop/Person.h"

#include <array>
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace stride {

/**
 * Venue attendance per person, venue type and day of week, as read from the
 * subpools community file. Build-time only: PopBuilder fills it, ContactDivider
 * adds the contact counts, and CopyToPools() moves what the simulation needs into
 * the pools (parallel to their member lists). It is released before the simulation
 * starts, so no per-person-per-day venue data stays resident (plan F12.3).
 */
class VenueAttendance
{
public:
        struct Record
        {
                unsigned int pool_id  = 0; ///< 0: no pool of this type on this day.
                unsigned int duration = 0;
                unsigned int contacts = 0;
        };

        /// One record per venue type (OtherHouse, RestoCafe, OtherPlace, Transport) per day.
        using Week = std::array<std::array<Record, 7>, 4>;

        /// Index of a venue type in Week.
        static std::size_t VenueIndex(ContactType::Id typ)
        {
                switch (typ) {
                case ContactType::Id::OtherHouse: return 0;
                case ContactType::Id::RestoCafe: return 1;
                case ContactType::Id::OtherPlace: return 2;
                case ContactType::Id::Transport: return 3;
                default:
                        throw std::runtime_error("VenueAttendance> not a venue type: " + ContactType::ToString(typ));
                }
        }

        /// The records of a person (created, all zero, on first access).
        Week& Ref(const Person& p) { return m_records[p.GetId()]; }

        /// The records of a person (all zero if the person has none).
        const Week& CRef(const Person& p) const
        {
                static const Week none{};
                const auto it = m_records.find(p.GetId());
                return it == m_records.end() ? none : it->second;
        }

        /// For every member of every venue pool, store the duration and contact count the
        /// member has on the pool's day. That is exactly what the transmission loop used to
        /// read per person, including for the rare pool whose members were read with
        /// different days (the pool then runs on the last day read).
        void CopyToPools(ContactPoolSys& poolSys) const
        {
                for (const auto typ : {ContactType::Id::OtherHouse, ContactType::Id::RestoCafe,
                                       ContactType::Id::OtherPlace, ContactType::Id::Transport}) {
                        const auto venue = VenueIndex(typ);
                        auto&      pools = poolSys.RefPools(typ);
                        for (size_t i = 1; i < pools.size(); i++) {
                                auto&      pool = pools[i];
                                const auto day  = pool.GetDayWeek();
                                for (size_t m = 0; m < pool.size(); m++) {
                                        Record r{};
                                        if (day < 7) {
                                                r = CRef(*pool[m])[venue][day];
                                        }
                                        pool.SetMemberAttendance(m, r.duration, r.contacts);
                                }
                        }
                }
        }

        /// Free the records once they have been copied to the pools.
        void Release() { std::unordered_map<unsigned int, Week>().swap(m_records); }

private:
        std::unordered_map<unsigned int, Week> m_records; ///< Keyed by person id.
};

} // namespace stride
