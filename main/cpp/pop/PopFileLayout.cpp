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
 *  Copyright 2026, Willem L
 */

/**
 * @file
 * Column layout of a population file: implementation.
 */

#include "PopFileLayout.h"

#include "pop/PopBuilder.h"
#include "util/StringUtils.h"

using namespace std;
using namespace stride::util;

namespace stride {

PopFileLayout::PopFileLayout() : m_separator(","), m_has_header(true), m_num_columns(0U) { m_column.fill(-1); }

PopFileLayout PopFileLayout::FromFirstLine(const string& firstLine)
{
        PopFileLayout layout;
        layout.m_separator = (firstLine.find(';') != string::npos) ? ";" : ",";

        const auto headers   = Split(firstLine, layout.m_separator);
        const auto name      = [&headers](size_t i) { return Trim(ToString(headers[i]), ToString('"')); };
        layout.m_num_columns = headers.size();

        // Positional layout, inferred from the third column name and the column count.
        const bool profession     = headers.size() > 2 && name(2) == "worker";
        const int  profession_adj = profession ? 2 : 0;
        if (profession) {
                layout.Set(Field::PersonId, 1);
                layout.Set(Field::Profession, 2);
        }
        layout.Set(Field::Age, 0);
        layout.Set(Field::Household, 1 + profession_adj);
        layout.Set(Field::School, 2 + profession_adj);
        layout.Set(Field::Workplace, 3 + profession_adj);
        layout.Set(Field::CommunityWeekend, 4 + profession_adj);
        layout.Set(Field::CommunityWeekday, 5 + profession_adj);

        // At most one extra column, and only when the column count matches exactly.
        if (headers.size() == static_cast<size_t>(7 + profession_adj)) {
                const auto extra_id = name(6 + profession_adj);
                if (extra_id == "household_cluster_id") {
                        layout.Set(Field::HouseholdCluster, 6 + profession_adj);
                } else if (extra_id == "collectivity_id") {
                        layout.Set(Field::Collectivity, 6 + profession_adj);
                }
        }
        return layout;
}

PopRecord PopFileLayout::Parse(const string& line, unsigned int defaultPersonId) const
{
        const auto values = Split(line, m_separator);
        const auto get    = [this, &values](Field f) -> unsigned int {
                const int c = m_column[static_cast<size_t>(f)];
                if (c < 0 || static_cast<size_t>(c) >= values.size()) {
                        return 0U;
                }
                return static_cast<unsigned int>(IntFromString(values[c]));
        };

        PopRecord r;
        r.age               = get(Field::Age);
        r.person_id         = Has(Field::PersonId) ? get(Field::PersonId) : defaultPersonId;
        r.profession        = get(Field::Profession);
        r.household         = get(Field::Household);
        r.school            = get(Field::School);
        r.workplace         = get(Field::Workplace);
        r.community_weekend = get(Field::CommunityWeekend);
        r.community_weekday = get(Field::CommunityWeekday);
        // The extra column is only read from rows with exactly as many columns as the header.
        if (values.size() == m_num_columns) {
                r.household_cluster = get(Field::HouseholdCluster);
                r.collectivity      = get(Field::Collectivity);
        }
        return r;
}

} // namespace stride
