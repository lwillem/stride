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
#include "healthcare/ConstantVaccine.h"
#include "util/RnMan.h"
#include "util/FileSys.h"
#include "util/LogUtils.h"
#include "util/StringUtils.h"

#include "util/Ptree.h"
#include <algorithm>
#include <numeric>
#include <utility>
#include <vector>

namespace stride {

using namespace stride::ContactType;
using namespace stride::util;
using namespace std;


ImmunitySeeder::ImmunitySeeder(const ptree& config, std::shared_ptr<util::RnMan> rnMan) : m_config(config), m_rn_man(rnMan) {}


void ImmunitySeeder::Seed(std::shared_ptr<Population> pop)
{
        // --------------------------------------------------------------
        // Population immunity (natural & vaccine induced immunity).
        // --------------------------------------------------------------
        const auto immunityProfile = m_config.get<std::string>("run.immunity_profile");
        Vaccinate("immunity", immunityProfile, pop->CRefPoolSys().CRefPools<Id::Household>(),pop);

        const auto vaccinationProfile = m_config.get<std::string>("run.vaccine_profile");
        if(vaccinationProfile == "Teachers"){
        	Vaccinate("vaccine", "Random", pop->CRefPoolSys().CRefPools<Id::School>(),pop);
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
                        // No clustering requested => use the O(N) path (same per-age quota, no
                        // rejection, no household-size bias; see RandomIndependent()).
                        if (linkProbability == 0) {
                                RandomIndependent(immunityPools, immunityDistribution, pop, false);
                        } else {
                                Random(immunityPools, immunityDistribution, linkProbability, pop, false);
                        }

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

			// linkProbability is 0 on this path: the Random/Cocoon profiles have no
			// clustering knob, so this always takes the fast path.
			if (linkProbability == 0) {
				RandomIndependent(immunityPools_selection, immunityDistribution, pop, true);
			} else {
				Random(immunityPools_selection, immunityDistribution, linkProbability, pop, true);
			}

		}
}



void ImmunitySeeder::RandomIndependent(const SegmentedVector<ContactPool>& pools,
                       vector<double>& immunityDistribution, std::shared_ptr<Population> pop,
                       const bool log_immunity)
{
        // No clustering is requested: bucket the candidates per age class, shuffle, and
        // take the first `quota`. Same exact per-age quota as Random(), O(N), no rejection.
        //
        // NOT equivalent in distribution to Random() at immunityLinkProbability == 0.
        // Random() draws a household uniformly and then one member, so it over-samples
        // people in small households (Dane WI, 50 %: 77 % immune when living alone, 23 % in
        // households of 8+). Here every candidate of an age class is equally likely.
        // Accepted as a correction, 2026-10-06 (refactoring plan, Phase 5c step 4).
        //
        // NOTE: deliberately NOT the household-pruning variant that was tried on
        // measles_usa_rm. That kept the with-replacement draw and added an isExhausted()
        // scan of the household per draw, which made the 99.91% case far slower rather
        // than faster (see doc/markdown/measles_usa_rm_merge_result.md section 3).

        const unsigned int maxAge = pop->GetMaxAge();
        auto&              logger = pop->RefEventLogger();

        // Keep the pool alongside each candidate so the [VACC] log lines stay identical
        // in content to the ones Random() emits.
        vector<vector<std::pair<Person*, const ContactPool*>>> byAge(maxAge + 1);
        for (const auto& c : pools) {
                for (const auto& p : c.GetPool()) {
                        if (!p->IsVaccinated()) {
                                byAge[p->GetAge()].emplace_back(p, &c);
                        }
                }
        }


        shared_ptr<ConstantVaccine::Properties> properties(
            new ConstantVaccine::Properties{"immunity", 1.0, 1.0, 1.0});

        for (unsigned int age = 0; age <= maxAge; age++) {
                auto& candidates = byAge[age];
                if (candidates.empty()) { continue; }

                // Identical to Random(): floor(unvaccinated count * rate) for this age class.
                const auto quota = static_cast<unsigned int>(
                    floor(static_cast<double>(candidates.size()) * immunityDistribution[age]));
                if (quota == 0) { continue; }

                // Sample without replacement by PARTIAL Fisher-Yates: only the first
                // `quota` positions are needed, so this is O(quota), not O(n log n).
                //
                // Deliberately NOT RnMan::Shuffle (std::shuffle over trng::lcg64): that
                // is unusably slow at this scale -- shuffling 6,971 elements did not
                // complete in 45 s, while the household-sized vectors Random() passes it
                // (2-3 elements) are unaffected, which is why the defect has stayed
                // hidden. SampleUniform01() is the draw used elsewhere in this file.
                auto&      rng = m_rn_man->at(0U);
                const auto n   = static_cast<unsigned int>(candidates.size());

                vector<unsigned int> order(n);
                iota(order.begin(), order.end(), 0U);

                for (unsigned int i = 0; i < quota; i++) {
                        unsigned int j = i + static_cast<unsigned int>(rng.SampleUniform01() * (n - i));
                        if (j >= n) { j = n - 1; } // guard against SampleUniform01() == 1.0
                        std::swap(order[i], order[j]);
                }

                for (unsigned int i = 0; i < quota; i++) {
                        Person*            p      = candidates[order[i]].first;
                        const ContactPool& p_pool = *candidates[order[i]].second;

                        auto vaccine = std::unique_ptr<Vaccine>(new ConstantVaccine(properties));
                        p->SetVaccine(vaccine);

                        if (log_immunity) {
                                logger->info("[VACC] {} {} {} {} {} {}", p->GetId(), p->GetAge(),
                                             ToString(p_pool.GetType()), p_pool.GetId(),
                                             p_pool.HasInfant(), 0);
                        }
                }
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
        auto       intGenerator       = m_rn_man->at(0U).GetUniformIntGenerator(0, poolsSize);
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
                m_rn_man->at(0U).Shuffle(indices);

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
                        if (m_rn_man->at(0).SampleUniform01() < (1 - immunityLinkProbability)) {
                                break;
                        }
                }
        }
}


} // namespace stride
