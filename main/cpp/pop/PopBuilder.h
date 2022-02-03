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
 *  Copyright 2017, 2018 Kuylen E, Willem L, Broeckhove J
 *  Copyright 2018, 2019 Jan Broeckhove and Bistromatics group.
 */

/**
 * @file
 * Initialize populations.
 */

#pragma once

#include <boost/property_tree/ptree_fwd.hpp>
#include <memory>
#include <spdlog/logger.h>

namespace stride {

class Population;
namespace util {
class RnMan;
}

/**
 * Initializes Population objects.
 */
class PopBuilder
{
public:
        /// Initializing constructor.
        /// \param config        Property_tree with general configuration settings.
        /// \param strideLogger  Logging.
        PopBuilder(const boost::property_tree::ptree& config,
                           std::shared_ptr<spdlog::logger> strideLogger = nullptr);

        /// Build Population and return it afterwards.
        /// The steps are:
        /// - Check input data.
        /// - Read persons from file and instantiate them.
        /// - Fill up the various type of contactpools.
        /// - Seed the population with contact survey participants.
        std::shared_ptr<Population> Build(std::shared_ptr<Population> pop);


private:
        /// Generates pop's individuals and return pop.
        std::shared_ptr<Population> MakePersons(std::shared_ptr<Population> pop);
        /// Optimized implementation of MakePersons
        std::shared_ptr<Population> MakePersonsOpt(std::shared_ptr<Population> pop);

        const boost::property_tree::ptree& m_config;        ///< Configuration property tree.
        std::shared_ptr<spdlog::logger>    m_stride_logger; /// Logger for build process.
};

// Transform strings to unsigned integers faster than FromString
inline unsigned long IntFromString(const std::string& s)
{
    return std::strtoul(s.c_str(), nullptr, 10);
}

} // namespace stride
