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
#include <numeric>
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

        // ============================== CHANGED START ==============================
        // Single on/off switch between the two complete workflows:
        //   - false (default): the ORIGINAL ImmunitySeeder workflow, unchanged -- both
        //     "immunity" and "vaccine" passes use the household-clustered Random(),
        //     no hesitancy filtering. An existing config with this key absent behaves
        //     exactly as before.
        //   - true: "mass immunize -> vaccinate iteratively" -- natural immunity is
        //     assigned in one independent per-age pass (RandomIndependent, no
        //     household structure), then children are vaccinated iteratively via the
        //     household-clustered Random(), against a hesitancy-filtered household
        //     subset. This is the Naive-Immunizer.R workflow.
        const auto massImmunize = m_config.get<bool>("run.mass_immunize", false);

        // Vaccine hesitancy: a subset of households refuse vaccination for everyone in
        // the family, regardless of any individual's immunity status. Read as a
        // flexible config parameter, same pattern as "*_link_probability", so its rate
        // can vary per simulation -- but only ever applied under the mass-immunize
        // workflow, and only for the "vaccine" pass against Household-type pools (the
        // "Teachers" path vaccinates via School pools, which this does not apply to).
        // A FIXED COUNT of households is excluded (round(rate * N), sampled without
        // replacement), not an independent per-household coin flip, matching the
        // sampling approach validated in R.
		SegmentedVector<ContactPool> eligibleImmunityPools(immunityPools);
        const bool isHouseholdPools = !immunityPools.empty() && immunityPools[0].GetType() == Id::Household;

        if (massImmunize && immunityType == "vaccine" && isHouseholdPools) {
                const auto hesitancyRate = m_config.get<double>("run." + ToLower(immunityType) + "_hesitancy_rate", 0.0);

                if (hesitancyRate > 0.0) {
                        const auto numHouseholds = static_cast<unsigned int>(immunityPools.size());
                        const auto numHesitant   = static_cast<unsigned int>(floor(hesitancyRate * numHouseholds + 0.5));
                        const auto numEligible   = numHouseholds - numHesitant;

                        vector<unsigned int> order(numHouseholds);
                        iota(order.begin(), order.end(), 0U);
                        m_rn_man->at(0U).Shuffle(order);

                        eligibleImmunityPools.clear();
                        for (unsigned int i = 0; i < numEligible; i++) {
                                eligibleImmunityPools.push_back(immunityPools[order[i]]);
                        }
                }
        }
        // =============================== CHANGED END ================================

        if (immunizationProfile == "AgeDependent") {
                        const auto   immunityFile = m_config.get<string>("run." + ToLower(immunityType) + "_distribution_file");
                        const ptree& immunity_pt  = FileSys::ReadPtreeFile(immunityFile);

                        for (unsigned int index_age = 0; index_age <= maxAge; index_age++) {
                                auto immunityRate = immunity_pt.get<double>("immunity.age" + std::to_string(index_age));
                                immunityDistribution.push_back(immunityRate);
                        }

                        // ============================== CHANGED START ==============================
                        // Under the mass-immunize workflow, natural immunity is assigned
                        // independently per age, ignoring household structure entirely. Under the
                        // original workflow (default), this is untouched -- Random() as always.
                        if (massImmunize && immunityType == "immunity") {
                                RandomIndependent(immunityDistribution, eligibleImmunityPools, pop);
                        } else {
                                linkProbability = m_config.get<double>("run." + ToLower(immunityType) + "_link_probability");
                                Random(eligibleImmunityPools, immunityDistribution, linkProbability, pop, false);
                        }
                        // =============================== CHANGED END ================================

		} else if(immunizationProfile == "Random" || immunizationProfile == "Cocoon") {

			// Initialize new ContactPool vector
			SegmentedVector<ContactPool> immunityPools_selection;

			// immunizationProfile == Random: copy all contact pools
			// immunizationProfile == Cocoon: copy all contact pools with an infant
			for (auto& c : eligibleImmunityPools) { // CHANGED: iterate the hesitancy-filtered set instead of immunityPools
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

        // ============================== CHANGED START ==============================
        // Same with-replacement household draw and same full-reshuffle-on-every-visit
        // mechanic as the ORIGINAL -- statistically this is the same sampling process,
        // not a restructured one. The only change: once a household is provably
        // exhausted (no remaining member is both unvaccinated and in an age bracket
        // with open quota), it's pruned from the draw pool so it can never be drawn
        // again. That's the only source of "wasted" draws the original had, and
        // pruning removes it without altering the sampling distribution.
        //
        // An earlier version of this fix visited each household once in a single
        // shuffled pass, fully resolving it before moving on. That was found, via R
        // validation against real population data, to systematically OVER-cluster
        // relative to the original at high target rates (household pair-correlation
        // excess +0.05 vs. the original's -0.02, consistent across 8 seeds) --
        // because forcing full resolution before moving on pushes households toward
        // all-or-nothing outcomes more than the original's interleaved redraws did.
        // This pruning approach was the one confirmed to match (-0.02 vs. -0.02).
        vector<unsigned int> activePools(pools.size());
        iota(activePools.begin(), activePools.end(), 0U);

        auto isExhausted = [&](const ContactPool& p_pool) {
                for (const auto& p : p_pool.GetPool()) {
                        if (!p->IsVaccinated() && populationBrackets[p->GetAge()] > 0) {
                                return false;
                        }
                }
                return true;
        };

        while (numImmune > 0 && !activePools.empty()) {
                const auto drawPos  = static_cast<unsigned int>(m_rn_man->at(0U).SampleUniform01() * activePools.size());
                const unsigned int poolIdx = activePools[drawPos];
                const ContactPool& p_pool  = pools[poolIdx];

                if (isExhausted(p_pool)) {
                        activePools[drawPos] = activePools.back();
                        activePools.pop_back();
                        continue;
                }

                // random pool, random order of members -- same as the original
                const auto           size = static_cast<unsigned int>(p_pool.GetPool().size());
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

                if (isExhausted(p_pool)) {
                        activePools[drawPos] = activePools.back();
                        activePools.pop_back();
                }
        }
        // =============================== CHANGED END ================================
}

// ============================== CHANGED START ==============================
void ImmunitySeeder::RandomIndependent(vector<double>& immunityDistribution,
                                       const SegmentedVector<ContactPool>& pools, std::shared_ptr<Population> pop)
{
		// retrieve the maximum age in the population
		unsigned int maxAge = pop->GetMaxAge();

		// Bucket every unvaccinated person by age, across the full population -- the
		// pools passed in are Households, and every person belongs to exactly one, so
		// their union is the whole population. No household structure or clustering is
		// used here; this mirrors the "adults are fixed immune/susceptible independent
		// of clustering" simplification.
		vector<vector<Person*>> byAge(maxAge + 1);
		for (auto& c : pools) {
				for (const auto& p : c.GetPool()) {
						if (!p->IsVaccinated()) {
								byAge[p->GetAge()].push_back(p);
						}
				}
		}

		shared_ptr<ConstantVaccine::Properties> properties(new ConstantVaccine::Properties{"immunity", 1.0,1.0,1.0});

		for (unsigned int age = 0; age <= maxAge; age++) {
				auto& candidates = byAge[age];
				if (candidates.empty()) continue;

				const auto quota = static_cast<unsigned int>(floor(candidates.size() * immunityDistribution[age]));
				if (quota == 0) continue;

				// Sample without replacement: shuffle indices, take the first `quota`.
				vector<unsigned int> order(candidates.size());
				iota(order.begin(), order.end(), 0U);
				m_rn_man->at(0U).Shuffle(order);

				for (unsigned int i = 0; i < quota; i++) {
						auto vaccine = std::unique_ptr<Vaccine>(new ConstantVaccine(properties));
						candidates[order[i]]->SetVaccine(vaccine);
				}
		}
}
// =============================== CHANGED END ================================


} // namespace stride