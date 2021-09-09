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
 *  Copyright 2020, Kuylen E, Willem L, Broeckhove J
 */

/**
 * @file
 * Implementation for the NonComplianceSeeder class.
 */

#include "ContactHeterogeneitySeeder.h"

#include "pop/Population.h"
#include "util/Exception.h"
#include "util/FileSys.h"
#include "util/RnMan.h"

#include <boost/property_tree/ptree.hpp>
#include <cassert>

using namespace boost::property_tree;
using namespace stride::ContactType;
using namespace stride::util;
using namespace std;

namespace stride {

ContactHeterogeneitySeeder::ContactHeterogeneitySeeder(const ptree& config, RnMan& rnMan) : m_config(config), m_rn_man(rnMan) {}

shared_ptr<Population> ContactHeterogeneitySeeder::Seed(shared_ptr<Population> pop)
{

	auto& population = *pop;

	boost::optional<string> community_contact_distribution = m_config.get_optional<string>("run.contact_distribution");

	if (community_contact_distribution) {

		// Get target overdispersion
		double community_contact_distribution_overdispersion = m_config.get<double>("run.contact_distribution_overdispersion");

		if (*community_contact_distribution == "Gamma") {
			// Use distribution with mean 1 and overdispersion = contact_distribution_overdispersion
			double shape = community_contact_distribution_overdispersion;
			double scale = 1 / shape;
			auto gamma_generator = m_rn_man.GetGammaGenerator(shape, scale, 0U);

			// Seed community contact factors
			for (size_t i = 0; i < population.size(); ++i) {
				population[i].SetIndividualContactFactor(gamma_generator());
			}

		}

	}


	// Seed non-compliance

	// Non-compliance by age

	// Retrieve the maximum age in the population
	unsigned int maxAge = pop->GetMaxAge();

	vector<double> nonComplianceDistribution;
	boost::optional<string> nonComplianceByAgeFile = m_config.get_optional<string>("run.non_compliance_by_age_file");
	if (nonComplianceByAgeFile) {
		const ptree& nonComplianceByAge_pt = FileSys::ReadPtreeFile(*nonComplianceByAgeFile);
		for (unsigned int index_age = 0; index_age <= maxAge; index_age++) {
			auto fractionNonCompliant = nonComplianceByAge_pt.get<double>("fraction_non_compliers.age" + std::to_string(index_age));
			nonComplianceDistribution.push_back(fractionNonCompliant);
		}
	} else {
		// Age does not determine non-compliance: set probability of being a non-complier to 1 for all ages
		for (unsigned int index_age = 0; index_age <= maxAge; index_age++) {
			nonComplianceDistribution.push_back(1);
		}
	}

	// Non-compliance in pools
	boost::optional<string> nonCompliancePooltype = m_config.get_optional<string>("run.non_compliance_pooltype");
	if (nonCompliancePooltype) {
		string nonComplianceType = m_config.get<string>("run.non_compliance_type");
		if (*nonCompliancePooltype == "Workplace") {
			SeedPools<Id::Workplace>(pop, nonComplianceType, nonComplianceDistribution);
		} else if (*nonCompliancePooltype == "Community") {
			SeedPools<Id::PrimaryCommunity>(pop, nonComplianceType, nonComplianceDistribution);
		}
	}

	return pop;
}

template <Id pooltype>
void ContactHeterogeneitySeeder::SeedPools(std::shared_ptr<Population> pop, string nonComplianceType, std::vector<double> nonComplianceByAge) {
	auto generator0to1 = m_rn_man.GetUniform01Generator(0U);

	if (nonComplianceType == "Random") {
		// Distribute a given number of non-compliers at random in the population

		// TODO : is there a way to check numPersonsInPools > targetNumNonCompliers?
		unsigned int targetNumNonCompliers = m_config.get<int>("run.num_non_compliers");
		unsigned int numNonCompliers = 0U;

		// Get all pools of given pooltype
		const SegmentedVector<ContactPool>& pools = pop->CRefPoolSys().CRefPools(pooltype);
		auto generatorPool = m_rn_man.GetUniformIntGenerator(0, static_cast<int>(pools.size()), 0U);

		while (numNonCompliers < targetNumNonCompliers) {
			// Draw a random pool
			const ContactPool& pool = pools[generatorPool()];
			const auto poolSize = static_cast<unsigned int>(pool.GetPool().size());
			if (poolSize <= 0) {
				continue;
			}
			// Draw a random person in the pool
			auto generatorPerson = m_rn_man.GetUniformIntGenerator(0, poolSize, 0U);
			Person& p = *pool[generatorPerson()];
			// Check if person has not been drawn before
			if (p.IsNonComplier(pooltype)) {
				continue;
			}
			// Check age distribution
			if (generator0to1() < nonComplianceByAge[p.GetAge()]) {
				// Set person to be a non-complier
				RegisterNonComplier(pop, p, pooltype);
				// Update number of remaining non-compliers to be found
				numNonCompliers++;
			}

		}
	} else if (nonComplianceType == "Hotspots") {
		// Non-compliance is dependant on the geographical location of one's household

		const auto hotspotsFile = m_config.get<string>("run.non_compliance_hotspots_file");
		const ptree& hotspots_pt = FileSys::ReadPtreeFile(hotspotsFile);

		for (auto district : hotspots_pt.get_child("hotspots")) {
			// Get total number of persons in hotspot
			unsigned int totalDistrictPop = 0U;
			// Save IDs of households in hotspot district
			vector<unsigned int> householdsInDistrict;
			for (auto household : district.second.get_child("households")) {
				unsigned int hhID = household.second.get_value<unsigned int>();
				householdsInDistrict.push_back(hhID);
				totalDistrictPop += static_cast<int>(pop->CRefPoolSys().CRefPools(Id::Household)[hhID].size());
			}
			// Calculate the number of non-compliers per 'hotspot' district
			// = number of individuals in households belonging to district * fraction non-compliers in district
			double fractionNonCompliers = district.second.get<double>("fraction_non_compliers");
			unsigned int targetNumNonCompliers = floor(totalDistrictPop * fractionNonCompliers);

			// Sample non-compliers in this district
			unsigned int numNonCompliers = 0U;

			auto generatorHH = m_rn_man.GetUniformIntGenerator(0, static_cast<int>(householdsInDistrict.size()), 0U);
			while (numNonCompliers < targetNumNonCompliers) {
				// Choose a random household in the hotspot district
				unsigned int hhID = householdsInDistrict[generatorHH()];
				const ContactPool& pool = pop->CRefPoolSys().CRefPools(Id::Household)[hhID];

				const auto poolSize = static_cast<unsigned int>(pool.GetPool().size());
				if (poolSize <= 0) {
					continue;
				}
				// Choose a random person in this household
				auto generatorPerson = m_rn_man.GetUniformIntGenerator(0, poolSize, 0U);
				Person& p = *pool[generatorPerson()];
				// Check if person has not been drawn before
				if (p.IsNonComplier(pooltype)) {
					continue;
				}
				// Check age distribution
				if (generator0to1() < nonComplianceByAge[p.GetAge()]) {
					// Set person to be a non-complier
					RegisterNonComplier(pop, p, pooltype);
					// Update the number of remaining non-compliers to be found
					numNonCompliers++;
				}

			}

		}

	} else if (nonComplianceType == "Pools") {
		// Make all individuals belonging to certain contact pools non-compliers

		const auto poolsFile = m_config.get<string>("run.non_compliance_pools_file");
		const ptree& pools_pt = FileSys::ReadPtreeFile(poolsFile);

		// Get IDs of pools in which members are non-compliant
		vector<int> nonCompliantPools;
		for (auto poolID: pools_pt.get_child("non_compliant_pools.pools")) {
			//unsigned int poolID = pool.second.get_value<unsigned int>();
			nonCompliantPools.push_back(poolID.second.get_value<unsigned int>());
		}

		// Set individuals belonging to pools to be non-compliant
		// FIXME pooltype for clustering may be different from pooltype for non-compliant actions!
		for (auto pool: pop->CRefPoolSys().CRefPools(Id::Household)) {
			const auto poolID = pool.GetId();
			// Check if pool is in nonCompliantPools
			if (find(nonCompliantPools.begin(), nonCompliantPools.end(), poolID) != nonCompliantPools.end()) {
				// Make each member a non-complier
				for (unsigned int i_p = 0; i_p < pool.size(); i_p++) {
					Person& p = *pool[i_p];
					RegisterNonComplier(pop, p, pooltype);
				}
			}
		}
	}
}


bool ContactHeterogeneitySeeder::RegisterNonComplier(std::shared_ptr<Population> pop, Person& p, Id pooltype)
{
	Population& population  = *pop;
	auto&       logger      = population.RefEventLogger();

	// Set person to be non-complier for the correct pooltypes
	p.SetNonComplier(pooltype);

	// Primary and Secondary community are treated as one when it comes to non-compliance
	if (pooltype == Id::PrimaryCommunity) {
		p.SetNonComplier(Id::SecondaryCommunity);
	}

	// Log person details
	logger->info("[NCOM] {} {} {} {} {} {} {} {} {} {} {} {} {} {}", p.GetId(), p.GetAge(),
									p.GetPoolId(Id::Household), p.IsNonComplier(Id::Household),
									p.GetPoolId(Id::K12School), p.IsNonComplier(Id::K12School),
									p.GetPoolId(Id::College), p.IsNonComplier(Id::College),
									p.GetPoolId(Id::Workplace), p.IsNonComplier(Id::Workplace),
									p.GetPoolId(Id::PrimaryCommunity), p.IsNonComplier(Id::PrimaryCommunity),
									p.GetPoolId(Id::SecondaryCommunity), p.IsNonComplier(Id::SecondaryCommunity));
	return true;

}

} // namespace stride
