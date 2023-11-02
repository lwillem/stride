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
 *  Copyright 2017, Kuylen E, Willem L, Broeckhove J
 */

/**
 * @file
 * Implementation for the Immunizer class.
 */

#include "ImmunitySeeder.h"

#include "pop/Person.h"
#include "pop/ConstantVaccine.h"
#include "util/RnMan.h"
#include "util/FileSys.h"
#include "util/LogUtils.h"
#include "util/StringUtils.h"

#include <boost/property_tree/ptree.hpp>
#include <numeric>
#include <vector>

namespace stride {

using namespace boost::property_tree;
using namespace stride::ContactType;
using namespace stride::util;
using namespace std;


ImmunitySeeder::ImmunitySeeder(const ptree& config, RnMan& rnMan) : m_config(config), m_rn_man(rnMan) {}


void ImmunitySeeder::Seed(std::shared_ptr<Population> pop)
{
        // --------------------------------------------------------------
        // Population immunity (natural & vaccine induced immunity).
        // --------------------------------------------------------------
        const auto immunityProfile = m_config.get<std::string>("run.immunity_profile");
        Vaccinate("immunity", immunityProfile, pop->CRefPoolSys().CRefPools<Id::Household>(),pop);

        const auto vaccinationProfile = m_config.get<std::string>("run.vaccine_profile");
        if(vaccinationProfile == "Teachers"){
        	Vaccinate("vaccine", "Random", pop->CRefPoolSys().CRefPools<Id::K12School>(),pop);
        } else {
        	Vaccinate("vaccine", vaccinationProfile, pop->CRefPoolSys().CRefPools<Id::Household>(),pop);
        }
}

void ImmunitySeeder::Vaccinate(const std::string& immunityType, const std::string& immunizationProfile,
                              const SegmentedVector<ContactPool>& immunityPools,std::shared_ptr<Population> pop)
{
        std::vector<double> immunityDistribution;
        double              linkProbability = 0;

        // retrieve the maximum age in the population
        unsigned int maxAge = pop->GetMaxAge();

        if (immunizationProfile == "AgeDependent") {
                        const auto   immunityFile = m_config.get<string>("run." + ToLower(immunityType) + "_distribution_file");
                        const ptree& immunity_pt  = FileSys::ReadPtreeFile(immunityFile);

                        linkProbability = m_config.get<double>("run." + ToLower(immunityType) + "_link_probability");

                        for (unsigned int index_age = 0; index_age <= maxAge; index_age++) {
                                auto immunityRate = immunity_pt.get<double>("immunity.age" + std::to_string(index_age));
                                immunityDistribution.push_back(immunityRate);
                        }
                        Random(immunityPools, immunityDistribution, linkProbability, pop, false);

		} else if(immunizationProfile == "Random" || immunizationProfile == "Cocoon") {

			// Initialize new ContactPool vector
			SegmentedVector<ContactPool> immunityPools_selection;

			// immunizationProfile == Random: copy all contact pools
			// immunizationProfile == Cocoon: copy all contact pools with an infant
			for (auto& c : immunityPools) {
				if(immunizationProfile == "Random" || c.HasInfant()){
					immunityPools_selection.push_back(c);
				}
			}

			// get immunity rate and
			const auto immunityRate     = m_config.get<double>("run." + ToLower(immunityType) + "_rate");
			const auto immunity_min_age = m_config.get<double>("run." + ToLower(immunityType) + "_min_age",0);
			const auto immunity_max_age = m_config.get<double>("run." + ToLower(immunityType) + "_max_age",maxAge);

			// Initialize a vector to store the immunity rate per age class [0-maxAge].
			for (unsigned int index_age = 0; index_age <= maxAge; index_age++) {
					if(index_age >= immunity_min_age && index_age <= immunity_max_age){
						immunityDistribution.push_back(immunityRate);
					} else{
						immunityDistribution.push_back(0);
					}
			}

			Random(immunityPools_selection, immunityDistribution, linkProbability, pop, true);

		}
}



void ImmunitySeeder::Random(const SegmentedVector<ContactPool>& pools, vector<double>& immunityDistribution,
                       double immunityLinkProbability,std::shared_ptr<Population> pop, const bool log_immunity)
{

		// retrieve the maximum age in the population
		unsigned int maxAge = pop->GetMaxAge();

		// Initialize a vector to count the population per age class [0-100].
        vector<double> populationBrackets(maxAge+1, 0.0);

        // Sampler for int in [0, pools.size()) and for double in [0.0, 1.0).
        const auto poolsSize          = static_cast<int>(pools.size());
        auto       intGenerator       = m_rn_man.GetUniformIntGenerator(0, poolsSize, 0U);
        auto       uniform01Generator = m_rn_man.GetUniform01Generator(0U);
        auto&      logger             = pop->RefEventLogger();

        // Count unvaccinated individuals per age class
        for (auto& c : pools) {
                for (const auto& p : c.GetPool()) {
                        if (!p->IsVaccinated()) {
                                populationBrackets[p->GetAge()]++;
                        }
                }
        }


        // Calculate the number of "new immune" individuals per age class.
        unsigned int numImmune = 0;
        for (unsigned int age = 0; age <= maxAge; age++) {
                populationBrackets[age] = floor(populationBrackets[age] * immunityDistribution[age]);
                numImmune += static_cast<unsigned int>(populationBrackets[age]);

        }

        //Simple vaccine immunity
        shared_ptr<ConstantVaccine::Properties> properties(new ConstantVaccine::Properties{"immunity", 1.0,1.0,1.0});

        // Sample immune individuals, until all age-dependent quota are reached.
        while (numImmune > 0) {
                // random pool, random order of members
                const ContactPool&   p_pool = pools[intGenerator()];
                const auto           size   = static_cast<unsigned int>(p_pool.GetPool().size());
                vector<unsigned int> indices(size);
                iota(indices.begin(), indices.end(), 0U);
                m_rn_man.Shuffle(indices, 0U);

                // loop over members, in random order
                for (unsigned int i_p = 0; i_p < size && numImmune > 0; i_p++) {
                        Person& p = *p_pool[indices[i_p]];
                        // if p is susceptible and his/her age class has not reached the quota => make immune
                        if (!p.IsVaccinated() && populationBrackets[p.GetAge()] > 0) {
                                auto vaccine = std::unique_ptr<Vaccine>(new ConstantVaccine(properties));
                                p.SetVaccine(vaccine);

                                populationBrackets[p.GetAge()]--;
                                numImmune--;
                                // TODO: check log_level
                                if(log_immunity){
                                	logger->info("[VACC] {} {} {} {} {} {}",
                                				 p.GetId(), p.GetAge(),ToString(p_pool.GetType()), p_pool.GetId(), p_pool.HasInfant(),0);
                                }
                        }
                        // random draw to continue in this pool or to sample a new one
                        if (uniform01Generator() < (1 - immunityLinkProbability)) {
                                break;
                        }
                }
        }
}


} // namespace stride
