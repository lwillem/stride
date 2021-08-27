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
 *  Copyright 2017, 2018 Kuylen E, Willem L, Broeckhove J
 */

/**
 * @file
 * Implementation for HealthSeeder class.
 */

#include "HealthSeeder.h"

#include "Health.h"
#include "disease/TransmissionProfile.h"
#include "pop/Population.h"
#include "util/Assert.h"
#include "util/RnHandler.h"
#include <boost/property_tree/ptree.hpp>
#include <omp.h>


using namespace boost::property_tree;
using namespace stride::util;
using namespace std;

namespace stride {

HealthSeeder::HealthSeeder(const boost::property_tree::ptree& diseasePt)
    : m_start_symptomatic(), m_time_asymptomatic(), m_time_infectious(), m_time_symptomatic(), m_probability_symptomatic(),
	  m_sympt_cnt_reduction_work_school(), m_sympt_cnt_reduction_community()
{
        GetDistribution(m_start_symptomatic, diseasePt, "disease.start_symptomatic");
        GetDistribution(m_time_asymptomatic, diseasePt, "disease.time_asymptomatic");
        GetDistribution(m_time_infectious, diseasePt, "disease.time_infectious");
        GetDistribution(m_time_symptomatic, diseasePt, "disease.time_symptomatic");

        AssertThrow((abs(m_start_symptomatic.back() - 1.0) < 1.e-10), "Error in start_symptomatic", nullptr);
        AssertThrow((abs(m_time_asymptomatic.back() - 1.0) < 1.e-10), "Error in time_asymptomatic", nullptr);
        AssertThrow((abs(m_time_infectious.back() - 1.0) < 1.e-10), "Error in time_infectious", nullptr);
        AssertThrow((abs(m_time_symptomatic.back() - 1.0) < 1.e-10), "Error in time_symptomatic", nullptr);


        // load age-specific probability to be symptomatic
        unsigned int maxAge = 110; //TODO
		for (unsigned int index_age = 0; index_age <= maxAge; index_age++) {
				auto probabilitySymptomatic = diseasePt.get<double>("disease.prob_symptomatic.age" + std::to_string(index_age),1);
				m_probability_symptomatic.push_back(probabilitySymptomatic);
		}


        m_sympt_cnt_reduction_work_school = diseasePt.get<double>("disease.sympt_cnt_reduction_work_school",1.0);
        m_sympt_cnt_reduction_community   = diseasePt.get<double>("disease.sympt_cnt_reduction_community",1.0);

}

void HealthSeeder::GetDistribution(vector<double>& distribution, const ptree& rootPt, const string& xmlTag)
{
        const boost::property_tree::ptree& subtree = rootPt.get_child(xmlTag);
        for (const auto& tree : subtree) {
                distribution.push_back(tree.second.get<double>(""));
        }
}

unsigned short int HealthSeeder::Sample(const vector<double>& distribution, double random01)
{
        auto ret = static_cast<unsigned short int>(distribution.size());
        for (unsigned short int i = 0; i < distribution.size(); i++) {
                if (random01 <= distribution[i]) {
                        ret = i;
                        break;
                }
        }
        return ret;
}

void HealthSeeder::Seed(const std::shared_ptr<stride::Population>& pop, const TransmissionProfile& transProfile, vector<util::RnHandler>& handlers)
{
        auto& population = *pop;

#pragma omp parallel num_threads(handlers.size())
        {
                auto& gen01 = handlers[static_cast<size_t>(omp_get_thread_num())];
#pragma omp for
                for (size_t i = 0; i < population.size(); ++i) {

                		// initiate start for symptomatic and infectious period
                		auto startSymptomatic          = 0;
                        auto startInfectiousness       = 0;

                        // sample from given distribution, but limit "start infectiousness" to day 1 (= one day after infection)
						while(startInfectiousness < 1){
							startSymptomatic          = Sample(m_start_symptomatic, gen01());
							startInfectiousness       = startSymptomatic - Sample(m_time_asymptomatic, gen01());
						}

                        const auto timeInfectious      = Sample(m_time_infectious, gen01());
                        auto timeSymptomatic           = Sample(m_time_symptomatic, gen01());


                        const bool isSymptomatic = gen01() <= m_probability_symptomatic[population[i].GetAge()];
                        if(!isSymptomatic){
                        	timeSymptomatic = 0;
                        }

						double relative_susceptibility = transProfile.GetIndividualSusceptibility(population[i].GetAge());

                        population[i].GetHealth() =
                            Health(startInfectiousness, startSymptomatic, timeInfectious, timeSymptomatic,
                            		m_sympt_cnt_reduction_work_school,m_sympt_cnt_reduction_community,
                            		relative_susceptibility);
                }
        }
}

} // namespace stride
