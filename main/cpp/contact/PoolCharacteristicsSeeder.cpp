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
 *  Copyright 2023
 */

/**
 * @file
 * Implementation for the Ventilation Heterogeneity class and Ventilation non compliance.
 */

#include "PoolCharacteristicsSeeder.h"

#include "contact/ContactType.h"
#include "pop/Population.h"
#include "contact/ContactPoolSys.h"
#include "util/Exception.h"
#include "util/FileSys.h"
#include "util/RnMan.h"

#include "util/Ptree.h"
#include <cassert>
#include <cmath>

using namespace stride::util;
using namespace std;

namespace stride {

using namespace ContactType;

PoolCharacteristicsSeeder::PoolCharacteristicsSeeder(const ptree& config, std::shared_ptr<util::RnMan> rnMan) : m_config(config), m_rn_man(rnMan) {}

shared_ptr<Population> PoolCharacteristicsSeeder::Seed(shared_ptr<Population> pop, const ptree& poolCharacteristicsPt)
{
	
	auto uniform01Number= m_rn_man->at(0U).SampleUniform01();

	auto& population = *pop;

	auto& poolSys = population.RefPoolSys();

	std::cout << "Start PoolCharacteristicsSeeder." << std::endl;


	// Ventilation
	std::optional<string> ventilation_distribution = m_config.get_optional<string>("run.ventilation_distribution");
	double ventilation_distribution_overdispersion;
	// Get target overdispersion
	if (ventilation_distribution){
		ventilation_distribution_overdispersion = m_config.get<double>("run.ventilation_distribution_overdispersion");
	}

	for (ContactType::Id typ : ContactType::IdList) {
		std::string typString = ToString(typ);
		double ventilationPoolTypeAverage = poolCharacteristicsPt.get<double>("pool_characteristics.ventilation_reduction." + typString,0.1);	
		if (ventilation_distribution) {
			if (*ventilation_distribution == "Gamma") {
				double shape = ventilation_distribution_overdispersion;
				double scale = ventilationPoolTypeAverage / shape;
				auto gamma_generator = m_rn_man->at(0U).GetGammaGenerator(shape, scale);
				for (auto& pool: poolSys.RefPools(typ)) {
					double pool_ventilation_probability = gamma_generator();
					pool.SetVentilation(pool_ventilation_probability);
				}
			}
		} else {
		for (auto& pool: poolSys.RefPools(typ)) {
			pool.SetVentilation(ventilationPoolTypeAverage);
		};
		}
	}

	// air_mass, pool_duration => needed for airborne transmision
	std::vector<double> average_area_per_person_vector;
    std::vector<double> variability_area_vector;
    std::vector<double> minimum_area_vector;
    std::vector<double> ceiling_height_vector;
	std::vector<int> pool_duration_vector;
	double average_area_per_person;
	double variability_area;
	double minimum_area;
	double ceiling_height;
	int pool_duration;
	for (ContactType::Id typ: ContactType::IdList) {
		if (typ != Id::Household && typ != Id::CommunityWeekend && typ != Id::CommunityWeekday && typ != Id::HouseholdCluster) {
			std::string typString = ToString(typ);
			if (typ == Id::School){
				unsigned int maxAge = population.GetMaxAge();
				for (unsigned int index_age = 0; index_age <= maxAge; index_age++) {
					double input_average_area_per_person = poolCharacteristicsPt.get<double>("pool_characteristics.average_area_per_person." + typString + ".age" + std::to_string(index_age), 1.0);
					double input_variability_area = poolCharacteristicsPt.get<double>("pool_characteristics.variability_area." + typString + ".age" + std::to_string(index_age), 1.0);
					double input_minimum_area = poolCharacteristicsPt.get<double>("pool_characteristics.minimum_area." + typString + ".age" + std::to_string(index_age), 70.0);
					int input_pool_duration = poolCharacteristicsPt.get<int>("pool_characteristics.pool_duration." + typString + ".age" + std::to_string(index_age), 280);
					average_area_per_person_vector.push_back(input_average_area_per_person);
					variability_area_vector.push_back(input_variability_area);
					minimum_area_vector.push_back(input_minimum_area);
					pool_duration_vector.push_back(input_pool_duration);
				}
			} else {
				average_area_per_person = poolCharacteristicsPt.get<double>("pool_characteristics.average_area_per_person." + typString,1.0);
				variability_area = poolCharacteristicsPt.get<double>("pool_characteristics.variability_area." + typString,1.0);
				minimum_area = poolCharacteristicsPt.get<double>("pool_characteristics.minimum_area." + typString,70.0);
				pool_duration = poolCharacteristicsPt.get<int>("pool_characteristics.pool_duration." + typString, 350);
			}

			if (typ == Id::Workplace){
				unsigned int maxProfession = 2;
				for (unsigned int index_profession = 0; index_profession <= maxProfession; index_profession++) {
					double input_ceiling_height = poolCharacteristicsPt.get<double>("pool_characteristics.ceiling_height." + typString + ".profession" + std::to_string(index_profession),3.0);
					ceiling_height_vector.push_back(input_ceiling_height);
				}
			} else {
				ceiling_height = poolCharacteristicsPt.get<double>("pool_characteristics.ceiling_height." + typString, 3.0);
			}
			
			for (size_t i = 1; i < poolSys.RefPools(typ).size(); i++) {
				auto& pool = poolSys.RefPools(typ)[i];
				const auto& pMembers = pool.m_members;
				const auto  pSize    = pMembers.size();
			
				if (typ == Id::School && pSize > 0){
					float age = pMembers[0]->GetAge();
					average_area_per_person = average_area_per_person_vector[age];
					variability_area = variability_area_vector[age];
					minimum_area = minimum_area_vector[age];
					pool_duration = pool_duration_vector[age];			
				}
				else if (typ == Id::Workplace && pSize > 0){
					unsigned int poolTypeSpecification =  pool[0]->GetProfession();
					ceiling_height = ceiling_height_vector[poolTypeSpecification];
				}
			
			double variatie = uniform01Number * variability_area;
    		double grootte = average_area_per_person + variatie;
    		grootte *= pSize;
    		grootte = max(grootte, minimum_area);
    		pool.SetAirMass(grootte*ceiling_height);
			if (typ == Id::School || typ == Id::Workplace) {	 
				for (size_t i_person1 = 0; i_person1 < pSize; i_person1++) {
					const auto p = pMembers[i_person1];
					for (int dayWeek = 1; dayWeek <= 5; ++dayWeek) {
							p->PoolDurations(typ)[dayWeek] = pool_duration;
					}
				}
			} else if (typ == Id::Collectivity){
				for (size_t i_person1 = 0; i_person1 < pSize; i_person1++) {
					const auto p = pMembers[i_person1];
					for (int dayWeek = 0; dayWeek <= 6; ++dayWeek) {
							p->PoolDurations(typ)[dayWeek] = pool_duration;
					}
				}
			}
				
	
		
			}
		}


        
	}
        
	return pop;

}

shared_ptr<Population> PoolCharacteristicsSeeder::NonCompliance(shared_ptr<Population> pop)
{

	auto& population = *pop;

	auto& logger = population.RefEventLogger();

	// Seed non-compliance

	// Non-compliance in pools
	std::optional<string> nonCompliancePooltype = m_config.get_optional<string>("run.non_compliance_pooltype");
	if (nonCompliancePooltype) {
		string nonComplianceType = m_config.get<string>("run.non_compliance_type");
		ContactType::Id nonComplianceTypeId = ToId(nonComplianceType);
		for (auto& pool:pop->RefPoolSys().RefPools(nonComplianceTypeId)){
			pool.SetNonComplier();
			logger->info("[CHG] {} ", pool.IsNonComplier());
		}
	}

	return pop;
}



} // namespace stride
