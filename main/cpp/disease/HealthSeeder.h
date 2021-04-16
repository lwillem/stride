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
 *  Copyright 2021, Willem L, Kuylen E, Libin P
 */

/**
 * @file
 * Header file for the HealthSeeder.
 */

#pragma once

#include <boost/property_tree/ptree_fwd.hpp>
#include <memory>
#include <string>
#include <vector>

#include "sim/HospitalisationConfig.h"

namespace stride {

class Population;
class TransmissionProfile;

namespace util {
class RnHandler;
}

/**
 * Seeds the population with Health data.
 */
class HealthSeeder
{
public:
        /// Constructor requires disease data.
        explicit HealthSeeder(const boost::property_tree::ptree& diseasePt);

        /// Seeds the population with Health data.
        void Seed(const std::shared_ptr<Population>& pop, const HospitalisationConfig &hc, const TransmissionProfile& transProfile, std::vector<util::RnHandler>& handlers);

private:
        /// Utility method to extract distribution from data in ptree.
        void GetDistribution(std::vector<double>& distribution, const boost::property_tree::ptree& rootPt,
                             const std::string& xmlTag);

        /// Sample for each of the health data item individually.
        unsigned short int Sample(const std::vector<double>& distribution, double random01);

private:
        std::vector<double> m_start_symptomatic;
        std::vector<double> m_time_asymptomatic;
        std::vector<double> m_time_infectious;
        std::vector<double> m_time_symptomatic;
        std::vector<double> m_probability_symptomatic;

        double             m_sympt_cnt_reduction_work_school;  ///< Proportional reduction of days in work/school pool when symptomatic
        double             m_sympt_cnt_reduction_community;    ///< Proportional reduction of days in the community pools when symptomatic

};

} // namespace stride
