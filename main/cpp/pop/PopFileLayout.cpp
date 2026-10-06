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

#include <cstdlib>
#include <map>
#include <stdexcept>

using namespace std;
using namespace stride::util;

namespace stride {

namespace {

using Field = PopFileLayout::Field;

/// Header name as compared: without surrounding quotes, blanks or a trailing CR, lower case.
string CleanName(const string& s) { return ToLower(Trim(Trim(s, " \t\r\n"), "\"")); }

/// Accepted header names, synonyms included; FieldName() gives the canonical one.
const map<string, Field>& NameTable()
{
        static const map<string, Field> table{{"age", Field::Age},
                                              {"person_id", Field::PersonId},
                                              {"worker", Field::Profession},
                                              {"profession", Field::Profession},
                                              {"household_id", Field::Household},
                                              {"school_id", Field::School},
                                              {"work_id", Field::Workplace},
                                              {"workplace_id", Field::Workplace},
                                              {"community_weekend", Field::CommunityWeekend},
                                              {"community_weekend_id", Field::CommunityWeekend},
                                              {"primary_community", Field::CommunityWeekend},
                                              {"community_weekday", Field::CommunityWeekday},
                                              {"community_weekday_id", Field::CommunityWeekday},
                                              {"secondary_community", Field::CommunityWeekday},
                                              {"household_cluster_id", Field::HouseholdCluster},
                                              {"collectivity_id", Field::Collectivity}};
        return table;
}

/// A value as it may appear in a data row: a number, empty, or NA.
bool IsDataValue(const string& token)
{
        const auto t = CleanName(token);
        if (t.empty() || t == "na") {
                return true;
        }
        char* end = nullptr;
        strtod(t.c_str(), &end);
        return end == t.c_str() + t.size();
}

} // namespace

const vector<Field>& PopFileLayout::RequiredFields()
{
        static const vector<Field> required{Field::Age,       Field::Household,        Field::School,
                                            Field::Workplace, Field::CommunityWeekend, Field::CommunityWeekday};
        return required;
}

string PopFileLayout::FieldName(Field f)
{
        switch (f) {
        case Field::Age: return "age";
        case Field::PersonId: return "person_id";
        case Field::Profession: return "worker";
        case Field::Household: return "household_id";
        case Field::School: return "school_id";
        case Field::Workplace: return "work_id";
        case Field::CommunityWeekend: return "community_weekend";
        case Field::CommunityWeekday: return "community_weekday";
        case Field::HouseholdCluster: return "household_cluster_id";
        case Field::Collectivity: return "collectivity_id";
        default: return "?";
        }
}

PopFileLayout::PopFileLayout() : m_separator(","), m_has_header(true) { m_column.fill(-1); }

PopFileLayout PopFileLayout::Positional(const vector<string>& names)
{
        // The layout used before columns were resolved by name: a "worker" third column means
        // person_id and worker come after age, and a seventh (or ninth) column is read only
        // when it is household_cluster_id or collectivity_id.
        PopFileLayout layout;
        const bool    profession     = names.size() > 2 && names[2] == "worker";
        const int     profession_adj = profession ? 2 : 0;
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
        if (names.size() == static_cast<size_t>(7 + profession_adj)) {
                const auto& extra_id = names[6 + profession_adj];
                if (extra_id == "household_cluster_id") {
                        layout.Set(Field::HouseholdCluster, 6 + profession_adj);
                } else if (extra_id == "collectivity_id") {
                        layout.Set(Field::Collectivity, 6 + profession_adj);
                }
        }
        return layout;
}

PopFileLayout PopFileLayout::FromFirstLine(const string& firstLine)
{
        const string separator = (firstLine.find(';') != string::npos) ? ";" : ",";
        const auto   tokens    = Split(firstLine, separator);

        vector<string> names;
        for (const auto& t : tokens) {
                names.push_back(CleanName(t));
        }

        // ------------------------------------------------------------------
        // No header: the first line is already a person. Read it positionally.
        // ------------------------------------------------------------------
        bool has_header = false;
        for (const auto& t : tokens) {
                has_header = has_header || !IsDataValue(t);
        }
        if (!has_header) {
                auto layout         = Positional(names);
                layout.m_separator  = separator;
                layout.m_has_header = false;
                for (size_t c = RequiredFields().size(); c < names.size(); ++c) {
                        layout.m_ignored.push_back("column " + to_string(c + 1));
                }
                return layout;
        }

        // ------------------------------------------------------------------
        // Header: resolve columns by name.
        // ------------------------------------------------------------------
        PopFileLayout layout;
        layout.m_separator = separator;
        const auto& table  = NameTable();
        vector<bool> recognised(names.size(), false);
        for (size_t c = 0; c < names.size(); ++c) {
                const auto it = table.find(names[c]);
                if (it == table.end()) {
                        continue;
                }
                if (layout.Has(it->second)) {
                        throw runtime_error("Population file header: columns '" +
                                            names[layout.Column(it->second)] + "' and '" + names[c] +
                                            "' both hold " + FieldName(it->second) + ".");
                }
                layout.Set(it->second, static_cast<int>(c));
                recognised[c] = true;
        }

        // ------------------------------------------------------------------
        // A field not found by name falls back to its positional column, if that
        // column's name was not recognised either.
        // ------------------------------------------------------------------
        const auto positional = Positional(names);
        for (size_t f = 0; f < static_cast<size_t>(Field::Count); ++f) {
                const auto field = static_cast<Field>(f);
                const int  c     = positional.Column(field);
                if (!layout.Has(field) && c >= 0 && static_cast<size_t>(c) < names.size() && !recognised[c]) {
                        layout.Set(field, c);
                        recognised[c] = true;
                        layout.m_positional.push_back("'" + names[c] + "' read as " + FieldName(field));
                }
        }
        for (size_t c = 0; c < names.size(); ++c) {
                if (!recognised[c]) {
                        layout.m_ignored.push_back("'" + names[c] + "'");
                }
        }

        for (const auto field : RequiredFields()) {
                if (!layout.Has(field)) {
                        throw runtime_error("Population file header has no column for " + FieldName(field) +
                                            ": " + firstLine);
                }
        }
        return layout;
}

PopRecord PopFileLayout::Parse(const string& line, unsigned int defaultPersonId) const
{
        const auto values = Split(line, m_separator);
        const auto get    = [this, &values](Field f) -> unsigned int {
                const int c = Column(f);
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
        r.household_cluster = get(Field::HouseholdCluster);
        r.collectivity      = get(Field::Collectivity);
        return r;
}

} // namespace stride
