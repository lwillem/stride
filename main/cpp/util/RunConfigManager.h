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
 * Produce run config ptree.
 */

#pragma once

#include "util/ConfigInfo.h"

#include "util/Ptree.h"
#include <string>
#include <vector>

namespace stride {
namespace util {

/**
 * Produce run config ptree.
 */
class RunConfigManager
{
public:
        /// Produce property tree for config with given name.
        static stride::util::ptree Create(const std::string& configName);

        /// Set of threadcounts to use for tests based an nomber of available OpenMP threads.
        static std::vector<unsigned int> CreateNumThreads(unsigned int maxNum = ConfigInfo::NumberAvailableThreads());

        /// Reconstitute property tree from string representation.
        static stride::util::ptree FromString(const std::string& s);

        /// Produce string representation of property tree.
        static std::string ToString(const stride::util::ptree& pt);

private:
        /// Produce Influenza config for scenario tests.
        static std::string CreateTestsInfluenza();

        /// Produce Measles config for scenario tests.
        static std::string CreateTestsMeasles();

        /// Produce Covid19 config for scenario tests.
        static std::string CreateTestsCovid19();
};

} // namespace util
} // namespace stride
