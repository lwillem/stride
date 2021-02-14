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
 *  Copyright 2020, Willem L
 */

/**
 * @file
 * Header for the PublicHealthAgency class.
 */

#pragma once

#include "contact/ContactPool.h"
#include "util/RnMan.h"
#include "util/RnHandler.h"
#include "util/FileSys.h"
#include "util/SegmentedVector.h"

#include <boost/property_tree/ptree_fwd.hpp>
#include <memory>

namespace stride {

class Calendar;
class Population;
class Sim;

/**
 * Sets intervention measures.
 */
class PublicHealthAgency
{
public:
	    /// Default constructor
	    PublicHealthAgency();

        /// Initializing PublicHealthAgency.
		void Initialize(const boost::property_tree::ptree& config);

        /// Public Health Strategy: look for contacts of infected cases and quarantine infected cases
		void PerformContactTracing(std::shared_ptr<Population> pop, util::RnHandler& rnHandler, const std::shared_ptr<Calendar> calendar);

		bool IsK12SchoolOff(unsigned int age, bool isPreSchoolOff, bool isPrimarySchoolOff, bool isSecondarySchoolOff, bool isCollegeOff);

		/// Is Contact tracing active today?
		bool IsContactTracingActive(const std::shared_ptr<Calendar> calendar) const;

        /// Trace one individual
        void Trace(Person& p_case,
                std::shared_ptr<Population> pop,
				util::RnHandler& rnHandler,
                const std::shared_ptr<Calendar> calendar);

private:
        //contact tracing configuration
        double m_detection_probability;   ///< Detection probability of symptomatic cases.
        double m_tracing_efficiency_household;  ///< Tracing probability for household members
        double m_tracing_efficiency_other;      ///< Tracing probability for non-household members
        unsigned int m_case_finding_capacity;  ///< Maximum number of symptomatic cases with contact tracing per day
        unsigned int m_delay_isolation_index;         ///< Number of days after symptom onset to perform a clinical test
        unsigned int m_delay_contact_tracing; ///< Number of days after clinical test to start contact tracing
		double m_test_false_negative;         ///< False negative rate of PCR tests
};

} // namespace stride
