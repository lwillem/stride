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

#include <boost/property_tree/ptree.hpp>
#include <cassert>

using namespace boost::property_tree;
using namespace stride::util;
using namespace std;

namespace stride {

using namespace ContactType;

PoolCharacteristicsSeeder::PoolCharacteristicsSeeder(const ptree& config, RnMan& rnMan) : m_config(config), m_rn_man(rnMan) {}

double berekenGrootte(int aantalPersonen, double gemiddeldeGrootte, double variabiliteit, double minimumGrootte) {
    // Bereken een willekeurige variatie rond het gemiddelde
    double variatie = ((rand() % 100) / 100.0) * variabiliteit;
    
    // Bereken de uiteindelijke grootte met variatie
    double grootte = gemiddeldeGrootte + variatie;
    
    // Vermenigvuldig de grootte met het aantal personen
    grootte *= aantalPersonen;

    // Controleer of de berekende grootte het minimum overschrijdt
    grootte = max(grootte, minimumGrootte);
    
    return grootte;
}

shared_ptr<Population> PoolCharacteristicsSeeder::Seed(shared_ptr<Population> pop)
{
	


	auto& population = *pop;

	auto& poolSys = population.RefPoolSys();

	auto& logger = population.RefEventLogger();

	// Functie om de grootte van de ruimte te berekenen op basis van het aantal personen

	std::cout << "Start PoolCharacteristicsSeeder." << std::endl;


	for (ContactType::Id typ: ContactType::IdList) {
		if (typ != Id::Household && typ != Id::PrimaryCommunity && typ != Id::SecondaryCommunity && typ != Id::HouseholdCluster) {
		
		for (size_t i = 1; i < poolSys.RefPools(typ).size(); i++) {
			auto& pool = poolSys.RefPools(typ)[i];
			const auto& pMembers = pool.m_members;
        	const auto  pSize    = pMembers.size();
			

			if (typ == Id::K12School){
				
    float age = pMembers[0]->GetAge();

			//	std::cout << "pool" << std::endl;
				if (age < 3) {
					double grootte = berekenGrootte(pool.m_members.size(), 5.0, 1.5, 20);
				//	std::cout << "grootte " << grootte << std::endl;
					pool.SetAirMass(grootte*3);
					 
						for (size_t i_person1 = 0; i_person1 < pSize; i_person1++) {
                			const auto p = pMembers[i_person1];
							for (int dayWeek = 1; dayWeek <= 5; ++dayWeek) {
       				 				p->PoolDurations(typ)[dayWeek] = 280;
							}
						}
				}
				
				 else if (age < 6) {
					double grootte = berekenGrootte(pool.m_members.size(), 3.75, 0.75, 45);
				//	std::cout << "grootte " << grootte << std::endl;
					pool.SetAirMass(grootte*3);
					
						for (size_t i_person1 = 0; i_person1 < pSize; i_person1++) {
                			const auto p = pMembers[i_person1];
							for (int dayWeek = 1; dayWeek <= 5; ++dayWeek) {
       				 				p->PoolDurations(typ)[dayWeek] = 280;
							}
						}
				}
				else if (age < 12) {
					double grootte = berekenGrootte(pool.m_members.size(), 2.5, 0.5, 45);
				//	std::cout << "grootte " << grootte << std::endl;
					pool.SetAirMass(grootte*3);
					
						for (size_t i_person1 = 0; i_person1 < pSize; i_person1++) {
                			const auto p = pMembers[i_person1];
							for (int dayWeek = 1; dayWeek <= 5; ++dayWeek) {
       				 				p->PoolDurations(typ)[dayWeek] = 280;
							}
						}
				}

				else {
					double grootte = berekenGrootte(pool.m_members.size(), 2, 2, 40);
				//	std::cout << "grootte " << grootte << std::endl;
					pool.SetAirMass(grootte*3);
					
						for (size_t i_person1 = 0; i_person1 < pSize; i_person1++) {
                			const auto p = pMembers[i_person1];
							for (int dayWeek = 1; dayWeek <= 5; ++dayWeek) {
       				 				p->PoolDurations(typ)[dayWeek] = 320;
							}
						}
				}
			}
			else if (typ == Id::College) {
					double grootte = berekenGrootte(pool.m_members.size(), 4, 2, 40);
				//	std::cout << "grootte " << grootte << std::endl;
					pool.SetAirMass(grootte*3);
					
						for (size_t i_person1 = 0; i_person1 < pSize; i_person1++) {
                			const auto p = pMembers[i_person1];
							for (int dayWeek = 1; dayWeek <= 5; ++dayWeek) {
       				 				p->PoolDurations(typ)[dayWeek] = 418;
							}
						}
			}
			
			else if (typ == Id::Workplace) {
				unsigned int poolTypeSpecification =  pool[0]->GetProfession();
				if (poolTypeSpecification == 1){
				double grootte = berekenGrootte(pool.m_members.size(), 1000, 400, 20);
				// std::cout << "grootte " << grootte << std::endl;
					pool.SetAirMass(grootte*10);
					
					 }
					else {
						double grootte = berekenGrootte(pool.m_members.size(), 7, 2, 15);
					//	std::cout << "grootte " << grootte << std::endl;
						pool.SetAirMass(grootte*3);
						
					}					
							for (size_t i_person1 = 0; i_person1 < pSize; i_person1++) {
                			const auto p = pMembers[i_person1];
							for (int dayWeek = 1; dayWeek <= 5; ++dayWeek) {
       				 				p->PoolDurations(typ)[dayWeek] = 440;
							}
						}

				}

			

			else if (typ == Id::OtherHouse) {	
						double grootte = berekenGrootte(pool.m_members.size(), 7, 10, 15);
					//	std::cout << "grootte " << grootte << std::endl;
						pool.SetAirMass(grootte*3);
						
							}
						
			else if (typ == Id::RestoCafe) {		
						double grootte = berekenGrootte(pool.m_members.size(), 2, 2, 70);
					//	std::cout << "grootte " << grootte << std::endl;
						pool.SetAirMass(grootte*3);
						
							}
						
			else if (typ == Id::Transport) {
						double grootte = berekenGrootte(pool.m_members.size(), 1, 0.5, 70);
					//	std::cout << "grootte " << grootte << std::endl;
						pool.SetAirMass(grootte*3);
						
							}
						
			else if (typ == Id::OtherPlace) {
						double grootte = berekenGrootte(pool.m_members.size(), 1, 0.5, 70);
					//	std::cout << "grootte " << grootte << std::endl;
						pool.SetAirMass(grootte*3);
						
							}		
			}
		}
// Ventilation


	boost::optional<string> ventilationFile = m_config.get_optional<string>("run.ventilation_file");    
    if (ventilationFile) {

        const auto fp = m_config.get<bool>("run.use_install_dirs") ? FileSys::GetDataDir() / (ventilationFile.value_or("")) : filesys::path(ventilationFile.value_or(""));
        
		ptree ventilationPt = FileSys::ReadPtreeFile(fp);
		
		for (ContactType::Id typ : ContactType::IdList) {
			std::string typString = ToString(typ);
			double ventilationInfo = ventilationPt.get<double>("ventilation_reduction." + typString,0);
			for (auto& pool: poolSys.RefPools(typ)) {
				pool.SetVentilation(ventilationInfo);
			};
		};


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
	boost::optional<string> nonCompliancePooltype = m_config.get_optional<string>("run.non_compliance_pooltype");
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
