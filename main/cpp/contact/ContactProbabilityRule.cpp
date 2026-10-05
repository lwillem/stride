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
 * Implementation for the ContactProbabilityRule.
 */

#include "ContactProbabilityRule.h"

#include "util/StringUtils.h"

#include <map>
#include <stdexcept>

namespace stride {
namespace ContactProbabilityRule {

using namespace std;

namespace {

const map<string, Id>& RuleMap()
{
        static const map<string, Id> rules{make_pair("MIN", Id::Min), make_pair("MEAN", Id::Mean)};
        return rules;
}

} // namespace

string ToString(Id r)
{
        static const map<Id, string> names{make_pair(Id::Min, "Min"), make_pair(Id::Mean, "Mean")};
        return names.at(r);
}

bool IsRule(const string& s) { return RuleMap().count(util::ToUpper(s)) == 1; }

Id ToRule(const string& s)
{
        const auto& rules = RuleMap();
        const auto  it    = rules.find(util::ToUpper(s));
        if (it == rules.end()) {
                throw runtime_error("ContactProbabilityRule::ToRule> unknown rule '" + s +
                                    "'; valid values are 'Min' and 'Mean'.");
        }
        return it->second;
}

} // namespace ContactProbabilityRule
} // namespace stride
