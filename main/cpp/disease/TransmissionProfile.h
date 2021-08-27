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
 *  Copyright 2019, Willem L, Kuylen E, Broeckhove J
 */

/**
 * @file
 * Header for the TransmissionProfile class.
 */

#pragma once

#include "pop/Person.h"
#include "util/RnHandler.h"

#include <boost/property_tree/ptree.hpp>
#include <vector>
#include <numeric>


namespace stride {

/**
 * Transmission probability from disease data.
 */
class TransmissionProfile
{
public:
	TransmissionProfile(): m_transmission_probability(0),
							m_transmission_probability_distribution("Constant"),
							m_transmission_probability_distribution_overdispersion(0),
							m_susceptibility_age(100),
							m_rel_transmission_asymptomatic(1),
							m_rel_susceptibility_children(1) {}

	/// Initialize.
	void Initialize(const boost::property_tree::ptree& configPT, const boost::property_tree::ptree& diseasePt);

	/// Return mean transmission probability.
	double GetHomogeneousProbability() const;

	/// Return mean age-specfici susceptibility adjustment factor.
	double GetSusceptibilityFactor() const;

	/// Return age-specific susceptibility adjustment factor.
	double GetIndividualSusceptibility(unsigned int age) const;

	/// Return age-, health-, and person-specific transmission probability.
	double GetProbability(Person* p_infected, Person* p_susceptible) const;

	/// Draw individual transmission probability from distribution.
	double GetIndividualInfectiousness(util::RnHandler& generator) const;

private:
	double 						m_transmission_probability;

	std::string 				m_transmission_probability_distribution;
	double 						m_transmission_probability_distribution_overdispersion;

	std::vector<double>			m_susceptibility_age;

    double            			m_rel_transmission_asymptomatic; ///< Relative reduction of transmission for asymptomatic cases
    double             			m_rel_susceptibility_children; ///< Relative reduction of susceptibility for children vs. adults

};

} // namespace stride

