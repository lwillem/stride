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
 *  Copyright 2020, Willem L, Kuylen E, Broeckhove J
 *
 */

/**
 * @file
 * Implementation of ContactLogMode.
 */

#include "EventLogMode.h"

#include "util/StringUtils.h"

#include <map>

namespace stride {
namespace EventLogMode {

using namespace std;

string ToString(Id l)
{
        static map<Id, string> names{make_pair(Id::None, "None"),
        							 	 make_pair(Id::Incidence, "Incidence"),
										 make_pair(Id::Transmissions, "Transmissions"),
										 make_pair(Id::Participants, "Participants"),
										 make_pair(Id::All, "All")};
        return names.at(l);
}

Id ToMode(const string& s)
{
        static map<string, Id> modes{make_pair("NONE", Id::None),
        	                         make_pair("INCIDENCE", Id::Incidence),
        					         make_pair("TRANSMISSIONS", Id::Transmissions),
									 make_pair("CONTACTTRACING", Id::Transmissions), // use Transmission logging by default
									 make_pair("PARTICIPANTS", Id::Participants),
                                     make_pair("ALL", Id::All)};
        const std::string      t = util::ToUpper(s);
        return modes.at(t);
}

} // namespace InfectorLogMode
} // namespace stride
