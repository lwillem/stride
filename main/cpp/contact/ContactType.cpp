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
 *  Copyright 2017, Kuylen E, Willem L, Broeckhove J
 *  Copyright 2018, Jan Broeckhove and Bistromatics group.
 *  Copyright 2019, Willem L, Kuylen E, Broeckhove J
 */

/**
 * @file
 * Implementation of ContactPoolType.
 */

#include "ContactType.h"

#include <map>
#include <cctype>     // voor std::isspace
#include <algorithm>  // voor std::transform en std::find_if

namespace stride {
namespace ContactType {

using namespace std;

namespace {

void to_upper(string& s)
{
        std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
}

} // namespace

bool IsId(const string& s)
{
        static map<string, Id> ids{
            make_pair("Household", Id::Household),
            make_pair("School", Id::School),
            // make_pair("College", Id::College),
            make_pair("Workplace", Id::Workplace),
            make_pair("Community_Weekend", Id::CommunityWeekend),
            make_pair("Community_Weekday", Id::CommunityWeekday),
			make_pair("HouseholdCluster", Id::HouseholdCluster),
			make_pair("Collectivity", Id::Collectivity),
            make_pair("OtherHouse", Id::OtherHouse),
            make_pair("RestoCafe", Id::RestoCafe),
            make_pair("OtherPlace", Id::OtherPlace),
            make_pair("Transport", Id::Transport)

        };
        string t{s};
        to_upper(t);
        return (ids.count(t) == 1);
}

Id ToId(const string& s)
{
        static map<string, Id> ids{
            make_pair("Household", Id::Household),
            make_pair("School", Id::School),
            make_pair("Workplace", Id::Workplace),
            make_pair("Community_Weekend", Id::CommunityWeekend),
            make_pair("Community_Weekday", Id::CommunityWeekday),
			make_pair("HouseholdCluster", Id::HouseholdCluster),
			make_pair("Collectivity", Id::Collectivity),
            make_pair("OtherHouse", Id::OtherHouse),
            make_pair("RestoCafe", Id::RestoCafe),
            make_pair("OtherPlace", Id::OtherPlace),
            make_pair("Transport", Id::Transport)

        };

    string t{s};
    //to_upper(t);
    // Verwijder aanhalingstekens rondom de invoerstring
    t.erase(std::remove_if(t.begin(), t.end(), [](char c) { return c == '"'; }), t.end());

    if (ids.count(t) == 1) {
        return ids[t];
    } else {
        cerr << "DEBUG: s = " << s << ", t = " << t << ", ids contents:" << endl;
        for (const auto& entry : ids) {
            cerr << entry.first << " -> " << static_cast<int>(entry.second) << " (Match: " << (entry.first == t) << ")" << endl;
        }
        throw runtime_error("ContactType::ToId> not available: " + s);
    }

       // return (ids.count(t) == 1) ? ids[t] : throw runtime_error("ContactType::ToId> not available:" + s + " " + t );
       
}

string ToString(Id c)
{
        static map<Id, string> names{
            make_pair(Id::Household, "Household"),
            make_pair(Id::School, "School"),
            make_pair(Id::Workplace, "Workplace"),
            make_pair(Id::CommunityWeekend, "CommunityWeekend"),
            make_pair(Id::CommunityWeekday, "CommunityWeekday"),
			make_pair(Id::HouseholdCluster, "HouseholdCluster"),
			make_pair(Id::Collectivity, "Collectivity"),
            make_pair(Id::OtherHouse, "OtherHouse"),
            make_pair(Id::RestoCafe, "RestoCafe"),
            make_pair(Id::OtherPlace, "OtherPlace"),
            make_pair(Id::Transport, "Transport")
        };
        return (names.count(c) == 1) ? names[c] : throw runtime_error("ContactType::ToString> not available:");
}

} // namespace ContactType
} // namespace stride
