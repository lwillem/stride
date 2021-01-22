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
 *  Copyright 2019, Kuylen E, Willem L, Broeckhove J
 */

/**
 * @file
 * Implementation for the TransmissionProfile class.
 */

#include "TransmissionProfile.h"

#include "util/StringUtils.h"

#include <cmath>
#include <boost/math/distributions/gamma.hpp>

namespace stride {

using namespace std;
using namespace boost::property_tree;
using namespace stride::util;

void TransmissionProfile::Initialize(const ptree& configPt, const ptree& diseasePt)
{
    // 1. setup general transmission aspects
    m_rel_transmission_asymptomatic   = diseasePt.get<double>("disease.rel_transmission_asymptomatic", 1);
    m_rel_susceptibility_children     = diseasePt.get<double>("disease.rel_susceptibility_children", 1);

    // 2. setup transmission probability: with a given R0 or a given mean transmission probability
    // Use boost:optional to check which parameters are available in the config file
    boost::optional<double> transmission_probability_as_input = configPt.get_optional<double>("run.transmission_probability");
    boost::optional<double> r0_as_input = configPt.get_optional<double>("run.r0");

    // If available, use mean transmission probability as input (dominates an input value for R0)
    if (transmission_probability_as_input) {
    		m_transmission_probability = *transmission_probability_as_input;

    	// Otherwise, calculate mean transmission probability from given R0
    } else {
    		double r0 = 0; // set R0 to 0 if missing
    	    	if (r0_as_input) {
    	    		r0 = *r0_as_input;
    	    	}

    	    	const auto b0 = diseasePt.get<double>("disease.transmission.b0");
    	    	const auto b1 = diseasePt.get<double>("disease.transmission.b1");
    	    	const auto b2 = diseasePt.get<double>("disease.transmission.b2");

    	    	// if linear model is fitted to simulation data:
    	    	if(b2 == 0){
    	    		// Expected(R0) = (b0 + b1*transm_prob).
    	    		m_transmission_probability = (r0 - b0) / b1;

    	    	} else{
    	    		// if quadratic model is fitted to simulation data:
    	    		// Expected(R0) = (b0 + b1*transm_prob + b2*transm_prob^2).
    	    		// Find root
    	    		const auto a = b2;
    	    		const auto b = b1;
    	    		const auto c = b0 - r0;

    	    		// To obtain a real values (instead of complex)
    	    		if (r0 < (-(b * b) / (4 * a))) {
    	    			const double determ = (b * b) - 4 * a * c;
    	    			m_transmission_probability = (-b + sqrt(determ)) / (2 * a);
    	    		} else {
    	    			throw runtime_error("TransmissionProfile::Initialize> Illegal input values.");
    	    		}
    	  	}
    }

    // Check if age-dependent susceptibility vector is available
    // Otherwise, susceptibility adjustment factor for all ages is 1.
    boost::optional<std::string> susceptibility_by_age_as_input = configPt.get_optional<std::string>("run.disease_susceptibility_age");
    boost::optional<std::string> susceptibility_agecat_as_input = configPt.get_optional<std::string>("run.disease_susceptibility_agecat");
    if (susceptibility_by_age_as_input && susceptibility_agecat_as_input) {

    		auto susceptibility_string = Split(*susceptibility_by_age_as_input, ",");
    		auto agecat_string = Split(*susceptibility_agecat_as_input, ",");

    		if(susceptibility_string.size() != agecat_string.size()){
    			throw runtime_error("TransmissionProfile::Initialize> Illegal input values for susceptibility_agecat: should be open ended with size susceptibility_age.size() ");
    		}

    		for (unsigned int index_agecat = 0; index_agecat < agecat_string.size(); index_agecat++) {
    			auto age_min = stod(agecat_string[index_agecat]);
    			auto age_max = (index_agecat == agecat_string.size()-1) ?
									m_susceptibility_age.size() :
									stod(agecat_string[index_agecat+1]) ;

 				for (unsigned int index_age = age_min; index_age < age_max; index_age++) {
						m_susceptibility_age[index_age] = stod(susceptibility_string[index_agecat]);
				}
    		}
    } else {
    		for (unsigned int index_age = 0; index_age < m_susceptibility_age.size(); index_age++) {
    			m_susceptibility_age[index_age] = 1;
    		}
    }

    // Check whether transmission probability follows a distribution (otherwise it remains constant / age)
    boost::optional<string> t_prob_distribution = configPt.get_optional<string>("run.transmission_probability_distribution");
    if (t_prob_distribution) {
    		m_transmission_probability_distribution = *t_prob_distribution;
    		// Get target overdispersion
    		m_transmission_probability_distribution_overdispersion = configPt.get<double>("run.transmission_probability_distribution_overdispersion");
    }


}

double TransmissionProfile::GetHomogeneousProbability() const {
	return m_transmission_probability;
}

double TransmissionProfile::GetSusceptibilityFactor() const {
	double susceptibility_mean = accumulate(m_susceptibility_age.begin(), m_susceptibility_age.end(), 0.0) / m_susceptibility_age.size();

	return susceptibility_mean;
}

double TransmissionProfile::GetIndividualSusceptibility(unsigned int age) const {
	if (age < m_susceptibility_age.size()) {
		return m_susceptibility_age[age];
	} else {
		return m_susceptibility_age[m_susceptibility_age.size() - 1];
	}
}

double TransmissionProfile::GetProbability(Person* p_infected, Person* p_susceptible) const {
	// Get individual transmission probability of infector
	double transmission_probability_infector = p_infected->GetHealth().GetRelativeInfectiousness();

	// Adjustment for asymptomatic cases
	double adjustment_asymptomatic = (p_infected->GetHealth().IsSymptomatic()) ? 1 : m_rel_transmission_asymptomatic;

	// Adjustment based on age-specific susceptibility of susceptible person

	// deprecated... this binary option will be removed in the future
	double adjustment_susceptible_child = (p_susceptible->GetAge() < 18) ? m_rel_susceptibility_children : 1;

	double adjustment_susceptible_age = p_susceptible->GetHealth().GetRelativeSusceptibility();

	return transmission_probability_infector * adjustment_asymptomatic * adjustment_susceptible_child * adjustment_susceptible_age;
}

double TransmissionProfile::GetIndividualInfectiousness(RnHandler& generator) const {

	// If mean transmission probability is 0, return 0.
	// FIXME Is this ok?
	if (m_transmission_probability == 0) {
		return m_transmission_probability;
	}

	if (m_transmission_probability_distribution == "Constant") {
		return m_transmission_probability;
	} else if (m_transmission_probability_distribution == "Gamma") {
		// Generate truncated (between 0 and 1) gamma distribution
		// Based on script https://rdrr.io/cran/RGeode/src/R/rgammatr.R
		double shape = m_transmission_probability_distribution_overdispersion;
		double scale = m_transmission_probability / shape;

		boost::math::gamma_distribution<double> gamma_dist = boost::math::gamma_distribution<double>(shape, scale);

		double cdf1 = cdf(gamma_dist, 0.0);
		double cdf2 = cdf(gamma_dist, 1.0);

		double individual_probability = quantile(gamma_dist, cdf1 + generator() * (cdf2 - cdf1));

		return individual_probability;

	} else {
		return m_transmission_probability;
	}

}

} // namespace stride
