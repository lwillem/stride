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
 *  Copyright 2018, Kuylen E, Willem L, Broeckhove J
 */

/**
 * @file
 * Header file for the NonComplianceSeeder class.
 */

#pragma once

#include <boost/property_tree/ptree_fwd.hpp>
#include <memory>

#include "contact/ContactType.h"
#include "contact/AgeContactProfiles.h"


namespace stride {

class Population;
class Sim;

namespace util {
class RnMan;
}

/**
 * Seed contact heterogeneity in the population,
 * either through non-compliers to social distancing measures
 * or through adding a distribution to contact rates in the community pools
 */
class ContactDivider
{
public:
	/// Initialize Seeder.
	/// \param config 		Configuration parameters.
	/// \param rnMan			Random number manager.
	ContactDivider(const boost::property_tree::ptree& config, util::RnMan& rnMan);

  
    /// \param pop               Population.
    /// \param ageContactProfiles Age Contact Profiles.
    std::shared_ptr<Population> Divide(std::shared_ptr<Population> pop, const AgeContactProfiles& ageContactProfiles);

private:
    const boost::property_tree::ptree& m_config; ///< Run config.
    util::RnMan&                       m_rn_man; ///< Random number manager.
    
};

} // namespace stride