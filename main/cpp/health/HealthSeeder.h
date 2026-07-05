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
 *  Copyright 2024
 */

/**
 * @file
 * Header file for the HealthSeeder.
 */

#pragma once

#include "pop/Age.h"

#include "util/Ptree.h"
#include <memory>
#include <string>
#include <vector>

namespace stride {

class Population;
class TransmissionProfile;

namespace util {
class RnMan;
}

/**
 * Seeds the population with Health data.
 */
class HealthSeeder
{
public:
        /// Constructor requires overall and disease-specifc data.
        explicit HealthSeeder(const stride::util::ptree& runPt,
        						const stride::util::ptree& diseasePt);

        /// Seeds the population with Health data.
       void Seed(const std::shared_ptr<Population>& pop, const TransmissionProfile& transProfile, std::shared_ptr<util::RnMan> rnMan);

private:
        /// Utility method to extract distribution from data in ptree.
        void GetDistribution(std::vector<double>& distribution, const stride::util::ptree& rootPt,
                             const std::string& xmlTag);

        /// Sample for each of the health data item individually.
        unsigned short int Sample(const std::vector<double>& distribution, double random01);

        /// Get the hospitalisation probability for an age.
		double GetHospitalProbability(const int age) const { return m_hospital_probabilities[EffectiveAge(age)]; }

		/// Get the hospitalisation delay for an age.
		double GetHospitalDelay(const int age) const { return m_hospital_delays[EffectiveAge(age)]; }

		/// Get the hospitalisation length of stay.
		double GetHospitalLengthOfStay() const { return m_hospital_length_of_stay; }

private:
        std::vector<double> m_start_symptomatic;
        std::vector<double> m_time_asymptomatic;
        std::vector<double> m_time_infectious;
        std::vector<double> m_time_symptomatic;
        std::vector<double> m_probability_symptomatic;

        double             m_sympt_cnt_reduction_workplace_school; ///< Proportional reduction of days in work/school pool when symptomatic
        double             m_sympt_cnt_reduction_community;        ///< Proportional reduction of days in the community pools when symptomatic

        std::array<double, MaximumAge() + 1> m_hospital_probabilities;  ///< Hospitalisation probabilities per age.
        std::array<double, MaximumAge() + 1> m_hospital_delays;        ///< Hospitalisation delays per age.
        unsigned short int m_hospital_length_of_stay;		            ///< Hospital length of stay.
};

} // namespace stride
