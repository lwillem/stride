#include "ContactDivider.h"
#include "contact/AgeContactProfile.h"
#include "contact/ContactType.h"
#include "pop/Population.h"
#include "contact/ContactPoolSys.h"
#include "util/Exception.h"
#include "util/FileSys.h"
#include "util/RnMan.h"

#include <boost/property_tree/ptree.hpp>
#include <cassert>
#include <random>
#include <vector>
#include <algorithm>
#include <iostream>
#include <functional>

using namespace boost::property_tree;
using namespace stride::util;
using namespace std;

namespace stride {

using namespace ContactType;

ContactDivider::ContactDivider(const ptree& config, RnMan& rnMan) : m_config(config), m_rn_man(rnMan) {}

shared_ptr<Population> ContactDivider::Divide(shared_ptr<Population> pop, const AgeContactProfiles& ageContactProfiles)
{

	std::cout << "Start ContactDivider" << std::endl;
	auto& population = *pop;

	auto& poolSys = population.CRefPoolSys();

	auto& logger = population.RefEventLogger();

	auto uniform01Generator = m_rn_man.GetUniform01Generator(0U);
    
	for (size_t i = 0; i < population.size(); ++i) {
		auto &p = population[i];
		unsigned int age = p.GetAge();

		for (size_t day = 0; day < 7; day++){
	   		const AgeContactProfile& profile = (day == 0 || day == 6) ?
                                       ageContactProfiles[Id::PrimaryCommunity] :
                                       ageContactProfiles[Id::SecondaryCommunity];

			double reference_num_contacts_p{profile[EffectiveAge(static_cast<unsigned int>(age))]};
			int rounded_reference_num_contacts_p = static_cast<int>(round(reference_num_contacts_p));

			unsigned int idOtherHouse = p.CPoolIds(Id::OtherHouse)[day];
			unsigned int idRestoCafe = p.CPoolIds(Id::RestoCafe)[day];
			unsigned int idOtherPlace = p.CPoolIds(Id::OtherPlace)[day];
			unsigned int idTransport = p.CPoolIds(Id::Transport)[day];
		
			unsigned int sizeOtherHouse = poolSys.CRefPools(Id::OtherHouse)[idOtherHouse].size();
			unsigned int sizeRestoCafe = poolSys.CRefPools(Id::RestoCafe)[idRestoCafe].size();
			unsigned int sizeOtherPlace = poolSys.CRefPools(Id::OtherPlace)[idOtherPlace].size();
			unsigned int sizeTransport = poolSys.CRefPools(Id::Transport)[idTransport].size();
						
			unsigned int durationOtherHouse = p.PoolDurations(Id::OtherHouse)[day];
			unsigned int durationRestoCafe = p.PoolDurations(Id::RestoCafe)[day];
			unsigned int durationOtherPlace = p.PoolDurations(Id::OtherPlace)[day];
			unsigned int durationTransport = p.PoolDurations(Id::Transport)[day];
			
			unsigned int totalDuration = durationOtherHouse + durationRestoCafe + durationOtherPlace + durationTransport;

        	double probabilityOtherHouse = static_cast<double>(durationOtherHouse) / totalDuration;
    		double probabilityRestoCafe = static_cast<double>(durationRestoCafe) / totalDuration;
    		double probabilityOtherPlace = static_cast<double>(durationOtherPlace) / totalDuration;
    		double probabilityTransport = static_cast<double>(durationTransport) / totalDuration;

			std::vector<double> probabilities = {probabilityOtherHouse,probabilityRestoCafe,probabilityOtherPlace,probabilityTransport};
    		std::vector<unsigned int> maxContactsPerLocation = {sizeOtherHouse - 1,sizeRestoCafe - 1, sizeOtherPlace - 1, sizeTransport -1};

			// Maak een vector met indices van 0 tot probabilities.size() - 1
    		std::vector<unsigned int> indices(probabilities.size());
    		std::iota(indices.begin(), indices.end(), 0);
			
			// Initialiseer resultaten
    		std::vector<unsigned int> results(probabilities.size(), 0);

    		for (unsigned int i = 0; i < rounded_reference_num_contacts_p; ++i) {
        		// Bereken aangepaste kansen op basis van reeds toegewezen contacten
        		std::vector<double> preservedProbabilities;
        		for (unsigned int index : indices) {
            		if (maxContactsPerLocation[index] > 0) {
                		preservedProbabilities.push_back(probabilities[index]);
            		}
        		}

				// Controleer of er nog beschikbare categorieën zijn
				if (!preservedProbabilities.empty()) {
				// Normaliseer de kansen
				double sum = std::accumulate(preservedProbabilities.begin(), preservedProbabilities.end(), 0.0);
				std::vector<double> normalizedProbabilities;
				std::vector<double> cumulativeProbabilities;
				double cumulative = 0.0;
				for (double prob : preservedProbabilities) {
    				double normalizedProb = prob / sum;
    				normalizedProbabilities.push_back(normalizedProb);
    				cumulative += normalizedProb;
    				cumulativeProbabilities.push_back(cumulative);
				}

        		unsigned int selectedCategory = 0;
				for (size_t i = 0; i < cumulativeProbabilities.size(); ++i) {
    				if (uniform01Generator() < cumulativeProbabilities[i]) {
        				selectedCategory = i;
        				break;
    				}
				}
        		
            	// Wijs een contact toe aan de geselecteerde categorie
            	results[selectedCategory]++;
            	maxContactsPerLocation[selectedCategory]--;
        		}
    		}
   		
			p.PoolContacts(Id::OtherHouse)[day] = results[0];
			p.PoolContacts(Id::RestoCafe)[day] = results[1];
			p.PoolContacts(Id::OtherPlace)[day] = results[2];
			p.PoolContacts(Id::Transport)[day] = results[3];
        } 

	}

	return pop;
}

} // namespace stride
