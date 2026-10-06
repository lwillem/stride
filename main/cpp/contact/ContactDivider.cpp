#include "ContactDivider.h"
#include "contact/AgeContactProfile.h"
#include "contact/ContactType.h"
#include "pop/Population.h"
#include "contact/ContactPoolSys.h"
#include "util/Exception.h"
#include "util/FileSys.h"
#include "util/RnMan.h"

#include <cassert>
#include <random>
#include <vector>
#include <algorithm>
#include <iostream>
#include <functional>

using namespace stride::util;
using namespace std;

namespace stride {

using namespace ContactType;

ContactDivider::ContactDivider(std::shared_ptr<util::RnMan> rnMan) : m_rn_man(rnMan) {}

shared_ptr<Population> ContactDivider::Divide(shared_ptr<Population> pop, const AgeContactProfiles& ageContactProfiles)
{

	std::cout << "Start ContactDivider" << std::endl;
	auto& population = *pop;

	auto& poolSys = population.CRefPoolSys();

	auto uniform01Number= m_rn_man->at(0U).SampleUniform01();
    
	for (size_t i = 0; i < population.size(); ++i) {
		auto &p = population[i];
		unsigned int age = p.GetAge();
		auto &week = population.RefVenueAttendance().Ref(p);

		for (size_t day = 0; day < 7; day++){
	   		const AgeContactProfile& profile = (day == 0 || day == 6) ?
                                       ageContactProfiles[Id::CommunityWeekend] :
                                       ageContactProfiles[Id::CommunityWeekday];

			double reference_num_contacts_p{profile[EffectiveAge(static_cast<unsigned int>(age))]};
			unsigned int rounded_reference_num_contacts_p = static_cast<unsigned int>(round(reference_num_contacts_p));

			unsigned int idOtherHouse = week[0][day].pool_id;
			unsigned int idRestoCafe = week[1][day].pool_id;
			unsigned int idOtherPlace = week[2][day].pool_id;
			unsigned int idTransport = week[3][day].pool_id;
		
			unsigned int sizeOtherHouse = poolSys.CRefPools(Id::OtherHouse)[idOtherHouse].size();
			unsigned int sizeRestoCafe = poolSys.CRefPools(Id::RestoCafe)[idRestoCafe].size();
			unsigned int sizeOtherPlace = poolSys.CRefPools(Id::OtherPlace)[idOtherPlace].size();
			unsigned int sizeTransport = poolSys.CRefPools(Id::Transport)[idTransport].size();
						
			unsigned int durationOtherHouse = week[0][day].duration;
			unsigned int durationRestoCafe = week[1][day].duration;
			unsigned int durationOtherPlace = week[2][day].duration;
			unsigned int durationTransport = week[3][day].duration;
			
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
    				if (uniform01Number < cumulativeProbabilities[i]) {
        				selectedCategory = i;
        				break;
    				}
				}
        		
            	// Wijs een contact toe aan de geselecteerde categorie
            	results[selectedCategory]++;
            	maxContactsPerLocation[selectedCategory]--;
        		}
    		}
   		
			week[0][day].contacts = results[0];
			week[1][day].contacts = results[1];
			week[2][day].contacts = results[2];
			week[3][day].contacts = results[3];
        } 

	}

	// move durations and contacts to the venue pools, parallel to their members
	population.RefVenueAttendance().CopyToPools(population.RefPoolSys());

	return pop;
}

} // namespace stride
