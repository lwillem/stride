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
 *  Copyright 2023
 */

/**
 * @file
 * Header file for the Ventilation Heterogeneity class and Ventilation non compliance.
 */

#pragma once

#include "util/Ptree.h"
#include <memory>

#include "contact/ContactType.h"


namespace stride {

class Population;

namespace util {
class RnMan;
}

/**
 * Seed ventilation heterogeneity in the pools,
 * either non-compliance
 */
class PoolCharacteristicsSeeder
{
public:
	/// Initialize Seeder.
	/// \param config 		Configuration parameters.
	/// \param rnMan			Random number manager.
	PoolCharacteristicsSeeder(const stride::util::ptree& config, std::shared_ptr<util::RnMan> rnMan);

    /// Fill extra characteristics in the contactPoolSys for airborne transmission.
    /// \param pop               Population.
    std::shared_ptr<Population> Seed(std::shared_ptr<Population> pop, const stride::util::ptree& poolCharacteristicsPt);

    /// \param pop               Population.
    std::shared_ptr<Population> NonCompliance(std::shared_ptr<Population> pop);


private:
    const stride::util::ptree&            m_config; ///< Run config.
    std::shared_ptr<util::RnMan>          m_rn_man; ///< Random number manager.
};

} // namespace stride
