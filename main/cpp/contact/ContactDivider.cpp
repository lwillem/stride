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

	for (size_t i = 0; i < population.size(); ++i) {
		auto &p population[i];

		unsigned int age = p.GetAge();

		for (size_t day = 0; day < 7; day++){
			if (day == 0 || day == 6) {
				const AgeContactProfile& profile = ageContactProfiles[Id::PrimaryCommunity];
			}
			else {
				const AgeContactProfile& profile = ageContactProfiles[Id::SecondaryCommunity];
			}

			double reference_num_contacts_p{profile[EffectiveAge(static_cast<unsigned int>(age))]};

			unsigned int idOtherHouse = p.CPoolIds(Id::OtherHouse)[day];
			unsigned int idRestoCafe = p.CPoolIds(Id::RestoCafe)[day];
			unsigned int idOtherPlace = p.CPoolIds(Id::OtherPlace)[day];
			unsigned int idTransport = p.CPoolIds(Id::Transport)[day];

			unsigned int sizeOtherHouse = poolSys(Id::OtherHouse)[idOtherHouse].size();
			unsigned int sizeRestoCafe = poolSys(Id::RestoCafe)[idRestoCafe].size();
			unsigned int sizeOtherPlace = poolSys(Id::OtherPlace)[idOtherPlace].size();
			unsigned int sizeTransport = poolSys(Id::Transport)[idTransport].size();
			
			unsigned int durationOtherHouse = p.CPoolDurations(Id::OtherHouse)[day];
			unsigned int durationRestoCafe = p.CPoolDurations(Id::RestoCafe)[day];
			unsigned int durationOtherPlace = p.CPoolDurations(Id::OtherPlace)[day];
			unsigned int durationTransport = p.CPoolDurations(Id::Transport)[day];

			unsigned int totalDuration = durationOtherHouse + durationRestoCafe + durationOtherPlace + durationTransport;

        	double probabilityOtherHouse = static_cast<double>(durationOtherHouse) / totalDuration;
    		double probabilityRestoCafe = static_cast<double>(durationRestoCafe) / totalDuration;
    		double probabilityOtherPlace = static_cast<double>(durationOtherPlace) / totalDuration;
    		double probabilityTransport = static_cast<double>(durationTransport) / totalDuration;

			std::vector<double> probabilities = {probabilityOtherHouse,probabilityRestoCafe,probabilityOtherPlace,probabilityTransport};
    		std::vector<unsigned int> maxContactsPerLocation = {sizeOtherHouse - 1,sizeRestoCafe - 1, sizeOtherPlace - 1, sizeTransport -1};

			std::random_device rd;
    		std::mt19937 gen(rd());
    
    		std::vector<unsigned int> result;
    		std::multinomial_distribution<unsigned int> distribution(reference_num_contacts_p, probabilities.begin(), probabilities.end());

			// Continue generating until a suitable distribution is found
    		bool validDistribution = false;
    		while (!validDistribution) {
        	result.clear();
        	result = distribution(gen);

        	// Check if the distribution satisfies the maximum conditions
        	validDistribution = std::all_of(result.begin(), result.end(), [&](unsigned int contacts) {
            return contacts <= maxContactsPerLocation[&contacts - &result[0]];
        	});
    		}

			p.PoolContacts(Id::OtherHouse)[day] = result[0];
			p.PoolContacts(Id::RestoCafe)[day] = results[1];
			p.PoolContacts(Id::OtherPlace)[day] = results[2];
			p.PoolContacts(Id::Transport)[day] = results[3];

        } 

	}
        
	return pop;

} // namespace stride
