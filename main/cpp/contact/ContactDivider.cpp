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
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

using namespace boost::property_tree;
using namespace stride::util;
using namespace std;

namespace stride {

using namespace ContactType;

ContactDivider::ContactDivider(const ptree& config, RnMan& rnMan) : m_config(config), m_rn_man(rnMan) {}

shared_ptr<Population> ContactDivider::Divide(shared_ptr<Population> pop, const AgeContactProfiles& ageContactProfiles)
{
	auto& population = *pop;

	auto& poolSys = population.CRefPoolSys();

	auto& logger = population.RefEventLogger();

	// Initialisatie van RNG
    const gsl_rng_type* rngType;
    gsl_rng* rng;

    gsl_rng_env_setup();
    rngType = gsl_rng_default;
    rng = gsl_rng_alloc(rngType);

	// Aantal categorieën (locaties)
    const size_t numCategories = 4;
    
	for (size_t i = 0; i < population.size(); ++i) {
		auto &p = population[i];

		unsigned int age = p.GetAge();

		for (size_t day = 0; day < 7; day++){
	   		const AgeContactProfile& profile = (day == 0 || day == 6) ?
                                       ageContactProfiles[Id::PrimaryCommunity] :
                                       ageContactProfiles[Id::SecondaryCommunity];

			double reference_num_contacts_p{profile[EffectiveAge(static_cast<unsigned int>(age))]};

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

			// Resultaten voor elke dag
    		std::vector<unsigned int> results(numCategories);
			bool validDistribution = false;

        	// Blijf proberen totdat een geldige verdeling is verkregen
        	while (!validDistribution) {
            // Simuleer multinomiale verdeling
            gsl_ran_multinomial(rng, numCategories, 1, probabilities.data(), results.data());

            // Controleer of de verdeling voldoet aan de maximale contacten per locatie
            validDistribution = true;
            for (size_t j = 0; j < numCategories; ++j) {
                if (results[j] > maxContactsPerLocation[j]) {
                    validDistribution = false;
                    break;
                }
           	 }
        	}
        
			p.PoolContacts(Id::OtherHouse)[day] = results[0];
			p.PoolContacts(Id::RestoCafe)[day] = results[1];
			p.PoolContacts(Id::OtherPlace)[day] = results[2];
			p.PoolContacts(Id::Transport)[day] = results[3];

        } 

	}

	// Vrijgeven van resources
    gsl_rng_free(rng);
        
	return pop;
}

} // namespace stride
