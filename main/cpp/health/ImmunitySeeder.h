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
 * Header for the ImmunitySeeder class.
 */

#pragma once

#include "contact/ContactPool.h"
#include "pop/Population.h"
#include "util/RnMan.h"
#include "util/SegmentedVector.h"

#include <boost/property_tree/ptree_fwd.hpp>
#include <memory>

namespace stride {

/**
 * Seeds population w.r.t immunity (natural immunity, vaccination, ...)
 */
class ImmunitySeeder
{
public:
		/// Initializing ImmunitySeeder.
        explicit ImmunitySeeder(const boost::property_tree::ptree& config, util::RnMan& rnMan);

        /// Build the simulator.
        void Seed(std::shared_ptr<Population> pop);

private:
        /// Seed for vaccination/natural immunity.
        void Vaccinate(const std::string& immunityType, const std::string& immunizationProfile,
                       const util::SegmentedVector<ContactPool>& immunityPools, std::shared_ptr<Population> pop);

        /// Random immunization.
		void Random(const util::SegmentedVector<ContactPool>& pools, std::vector<double>& immunityDistribution,
					double immunityLinkProbability, std::shared_ptr<Population> pop, const bool log_immunity);


private:
		const boost::property_tree::ptree& m_config; ///< Run config.
		util::RnMan& m_rn_man; ///< Random number manager.
};

} // namespace stride
