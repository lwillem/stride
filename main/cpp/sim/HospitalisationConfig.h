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
 * Header for the HospitalisationConfig class.
 */

#pragma once

#include <vector>
#include <array>

#include <pop/Age.h>

namespace stride {

/**
 * Class to represent hospitalisation configuration (probability, delay).
 */
class HospitalisationConfig 
{
public:
        /// Default constructor.
        HospitalisationConfig();

        /// Constructor accepting age categories, and for each category the probability and delay to hospitalisation.
        HospitalisationConfig(std::vector<unsigned int> ageCategories, 
                std::vector<double> probabilities, std::vector<double> delays);

        /// Get the hospitalisation probability for an age.
        double GetProbability(const int age) const { return m_probabilities[EffectiveAge(age)-1]; }

        /// Get the hospitalisation delay for an age.
        double GetDelay(const int age) const { return m_delays[EffectiveAge(age)-1]; }

private:
        void NoHospitalisationInit();

private:
        std::array<double, MaximumAge()> m_probabilities; ///< Hospitalisation probabilities per age. 
        std::array<double, MaximumAge()> m_delays;        ///< Hospitalisation probabilities per age. 
};

} // namespace stride
