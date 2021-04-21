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
 * Implementation for the HospitalisationConfig class.
 */

#include "HospitalisationConfig.h"

namespace stride {

void HospitalisationConfig::NoHospitalisationInit()
{
    for (unsigned int i=0; i <= MaximumAge(); ++i) {
        m_probabilities[i] = 0.0;
        m_delays[i] = 0.0;
    }   
}

HospitalisationConfig::HospitalisationConfig() : m_probabilities(), m_delays() 
{
    NoHospitalisationInit();
}

HospitalisationConfig::HospitalisationConfig(
    std::vector<unsigned int> ageCategories,
    std::vector<double> probabilities, std::vector<double> delays)
    : m_probabilities(), m_delays()
{
    if (ageCategories.empty()) {
        NoHospitalisationInit();
    } else {
        if (ageCategories[0] > 0) {
            for (unsigned int i=0; i <= ageCategories[0]; ++i) {
                m_probabilities[i] = 0.0;
                m_delays[i] = 0.0;
            }
        }
        for (unsigned int i=0; i <= ageCategories.size() - 1; ++i) {
            unsigned int max = MaximumAge();
            if (i < ageCategories.size() - 1)
		max = ageCategories[i+1];
            for (unsigned int j=ageCategories[i]+1; j <= max; ++j) {
                m_probabilities[j] = probabilities[i];
                m_delays[j] = delays[i];
            } 
        }
    }
}

}
