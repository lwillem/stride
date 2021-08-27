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
 *  Copyright 2020, Libin P, Willem L
 */

/**
 * @file
 * Header for the UniversalTesting class.
 */

#pragma once

#include "PublicHealthAgency.h"
#include "PCRPool.h"
#include "contact/ContactPool.h"
#include "util/RnMan.h"
#include "util/RnHandler.h"
#include "util/SegmentedVector.h"
#include "util/FileSys.h"

#include <boost/property_tree/ptree_fwd.hpp>

#include <set>
#include <vector>

namespace stride {

class Calendar;
class Population;
class Sim;

class UniversalTesting 
{
public:
	    /// Default constructor
	    UniversalTesting();

		void Initialize(const boost::property_tree::ptree& config);

       /// Public Health Strategy: perform universal testing
       void PerformUniversalTesting(std::shared_ptr<Population> pop, util::RnHandler& rnHandler, const std::shared_ptr<Calendar> calendar, PublicHealthAgency& pha);

private:
        filesys::path m_unitest_planning_output_fn; ///> Filename to output the planning to
        //universal testing configuration
        std::string m_unitest_pool_allocation; ///< File that lists the pool to which households belong
        double m_unitest_fnr;             ///< False negative rate for pool testing (universal testing)
        unsigned int m_unitest_n_tests_per_day; ///< Number of PCR tests per day (universal testing)
        unsigned int m_unitest_pool_size; ///< Pool size (universal testing)
        double m_unitest_test_compliance; ///< Household compliance with testing (universal testing)
        double m_unitest_isolation_compliance; ///< Household compliance when isolated (universal testing)
        unsigned int m_unitest_isolation_delay; ///< Delay (in days) after which positive individuals are isolated (universal testing)
        unsigned int m_unitest_detectable_delay; ///< Delay (in days) after which positive individuals become PCR detectable (universal testing)
        std::string m_unitest_isolation_strategy; ///< Isolation strategy: isolate-pool/trace (universal testing)
        //universal testing planning
        std::vector<std::set<PCRPool>> m_unitest_planning; ///< Vector with at each element, a set of PCR pools, to be performed at one day
        unsigned int m_unitest_day_in_sweep; ///< The n-th day of the current universal testing sweep
};

} // namespace stride
