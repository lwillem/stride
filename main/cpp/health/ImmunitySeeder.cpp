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
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <string>
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
        bool                log_immunity    = false;

        // retrieve the maximum age in the population
        unsigned int maxAge = pop->GetMaxAge();

        // --------------------------------------------------------------
        // 1. Age-specific immunity/vaccination rates (no random numbers).
        // --------------------------------------------------------------
        if (immunizationProfile == "AgeDependent") {
                const auto   immunityFile = m_config.get<string>("run." + ToLower(immunityType) + "_distribution_file");
                const ptree& immunity_pt  = FileSys::ReadPtreeFile(immunityFile);

                linkProbability = m_config.get<double>("run." + ToLower(immunityType) + "_link_probability");

                for (unsigned int index_age = 0; index_age <= maxAge; index_age++) {
                        auto immunityRate = immunity_pt.get<double>("immunity.age" + std::to_string(index_age));
                        immunityDistribution.push_back(immunityRate);
                }

        } else if (immunizationProfile == "Random" || immunizationProfile == "Cocoon") {
                // linkProbability stays 0: the Random/Cocoon profiles have no clustering knob.
                log_immunity = true;

                const auto immunityRate     = m_config.get<double>("run." + ToLower(immunityType) + "_rate");
                const auto immunity_min_age = m_config.get<double>("run." + ToLower(immunityType) + "_min_age",0);
                const auto immunity_max_age = m_config.get<double>("run." + ToLower(immunityType) + "_max_age",maxAge);

                for (unsigned int index_age = 0; index_age <= maxAge; index_age++) {
                        if(index_age >= immunity_min_age && index_age <= immunity_max_age){
                                immunityDistribution.push_back(immunityRate);
                        } else{
                                immunityDistribution.push_back(0);
                        }
                }
        } else {
                // e.g. "None": nothing to do
                return;
        }

        // --------------------------------------------------------------
        // 2. Vaccine hesitancy (optional, "vaccine" pass over Household pools only).
        //    Splits the households into compliant and hesitant ones. Without
        //    run.vaccine_hesitancy_rate (or 0), no random numbers are drawn and all
        //    pools are used as before, so existing runs are unchanged.
        // --------------------------------------------------------------
        SegmentedVector<ContactPool> compliantPools;
        SegmentedVector<ContactPool> hesitantPools;
        SegmentedVector<ContactPool> hesitantClusterPools;
        const bool applyHesitancy = SelectHesitantHouseholds(immunityType, immunityPools, immunityDistribution, pop,
                                                             compliantPools, hesitantPools, hesitantClusterPools);

        // --------------------------------------------------------------
        // 3. Sample immune persons per age class.
        // --------------------------------------------------------------
        const auto immunize = [&](const SegmentedVector<ContactPool>& pools, std::vector<double>& distribution) {
                const SegmentedVector<ContactPool>* target = &pools;

                // immunizationProfile == Cocoon: only pools with an infant
                SegmentedVector<ContactPool> cocoonPools;
                if (immunizationProfile == "Cocoon") {
                        for (auto& c : pools) {
                                if (c.HasInfant()) {
                                        cocoonPools.push_back(c);
                                }
                        }
                        target = &cocoonPools;
                }

                // No clustering requested => use the O(N) path (same per-age quota, no
                // rejection, no household-size bias; see RandomIndependent()).
                if (linkProbability == 0) {
                        RandomIndependent(*target, distribution, pop, log_immunity);
                } else {
                        Random(*target, distribution, linkProbability, pop, log_immunity);
                }
        };

        if (!applyHesitancy) {
                immunize(immunityPools, immunityDistribution);
        } else {
                // compliant households: the rates as given
                immunize(compliantPools, immunityDistribution);

                // hesitant households: optionally a reduced uptake, relative to the rates above
                // (default 0 = refuse all). Clustered hesitant households can have their own
                // uptake (default: the same as the other hesitant households).
                const auto uptake        = m_config.get<double>("run.vaccine_hesitancy_uptake", 0.0);
                const auto uptakeCluster = m_config.get<double>("run.vaccine_hesitancy_cluster_uptake", uptake);
                for (const auto u : {uptake, uptakeCluster}) {
                        if (u < 0.0 || u > 1.0) {
                                throw runtime_error("ImmunitySeeder> vaccine_hesitancy_(cluster_)uptake must be in [0,1], got " +
                                                    to_string(u));
                        }
                }
                const auto immunizeReduced = [&](const SegmentedVector<ContactPool>& pools, double u) {
                        if (u > 0.0 && !pools.empty()) {
                                std::vector<double> reducedDistribution(immunityDistribution);
                                for (auto& rate : reducedDistribution) {
                                        rate *= u;
                                }
                                immunize(pools, reducedDistribution);
                        }
                };
                immunizeReduced(hesitantPools, uptake);
                immunizeReduced(hesitantClusterPools, uptakeCluster);
        }
}

bool ImmunitySeeder::SelectHesitantHouseholds(const std::string& immunityType,
                                              const SegmentedVector<ContactPool>& pools,
                                              const std::vector<double>& distribution,
                                              std::shared_ptr<Population> pop,
                                              SegmentedVector<ContactPool>& compliantPools,
                                              SegmentedVector<ContactPool>& hesitantPools,
                                              SegmentedVector<ContactPool>& hesitantClusterPools)
{
        // Only the vaccine pass is affected: natural immunity ("immunity" pass) is assigned
        // regardless of a household's attitude towards vaccination.
        if (immunityType != "vaccine") {
                return false;
        }

        const auto hesitancyRate = m_config.get<double>("run.vaccine_hesitancy_rate", 0.0);
        if (hesitancyRate <= 0.0) {
                return false;
        }
        if (hesitancyRate > 1.0) {
                throw runtime_error("ImmunitySeeder> run.vaccine_hesitancy_rate must be in [0,1], got " +
                                    to_string(hesitancyRate));
        }

        // Hesitancy is defined per household. The "Teachers" profile vaccinates via School
        // pools, where it does not apply. Pool 0 is the placeholder pool of the right type,
        // so the type of the container can be read from it.
        if (pools.empty() || pools[0].GetType() != Id::Household) {
                return false;
        }

        // ------------------------------------------------------------------
        // The hesitancy rate is the share of ALL (non-empty) households that is hesitant,
        // e.g. an adult vaccine-hesitancy estimate applied per household. Only households
        // with at least one member whose age has a vaccination rate > 0 (with the child
        // file: a member aged 0-18) are "eligible": hesitancy only changes anything there.
        // ------------------------------------------------------------------
        const auto isEligibleAge = [&](const Person* p) {
                const auto age = static_cast<unsigned int>(p->GetAge());
                return age < distribution.size() && distribution[age] > 0.0;
        };

        vector<bool> isEligible(pools.size(), false);
        unsigned int numEligible   = 0;
        unsigned int numHouseholds = 0;
        for (unsigned int i = 0; i < pools.size(); i++) {
                if (pools[i].size() > 0) {
                        numHouseholds++;
                }
                for (const auto& p : pools[i].GetPool()) {
                        if (isEligibleAge(p)) {
                                isEligible[i] = true;
                                numEligible++;
                                break;
                        }
                }
        }

        const auto numHesitant =
            static_cast<unsigned int>(floor(hesitancyRate * static_cast<double>(numHouseholds) + 0.5));

        // Sample `count` elements of `candidates` uniformly without replacement (partial
        // Fisher-Yates, O(count)) and mark them hesitant. Deliberately not RnMan::Shuffle,
        // which is far too slow at this scale (see RandomIndependent below).
        auto&        rng = m_rn_man->at(0U);
        vector<bool> isHesitant(pools.size(), false);
        const auto   markHesitant = [&](vector<unsigned int>& candidates, unsigned int count) {
                const auto n = static_cast<unsigned int>(candidates.size());
                count        = std::min(count, n);
                for (unsigned int i = 0; i < count; i++) {
                        unsigned int j = i + static_cast<unsigned int>(rng.SampleUniform01() * (n - i));
                        if (j >= n) { j = n - 1; } // guard against SampleUniform01() == 1.0
                        std::swap(candidates[i], candidates[j]);
                        isHesitant[candidates[i]] = true;
                }
                return count;
        };

        // ------------------------------------------------------------------
        // Optional clustering: concentrate hesitant households around a few contact pools
        // (e.g. schools), as in outbreaks centred on a handful of schools or congregations.
        //   run.vaccine_hesitancy_cluster_pooltype  pool type ("School", "Community_Weekend",
        //                                           ...); "None" or absent = no clustering
        //   run.vaccine_hesitancy_num_clusters      number of such pools, drawn at random
        //   run.vaccine_hesitancy_cluster_rate      share of the eligible households linked to
        //                                           those pools (any member attends) that is
        //                                           hesitant
        // The total number of hesitant households stays round(rate * all households): the
        // clustered ones are placed first, the remainder is spread at random over all other
        // households (if the clusters need more than the total, the total wins).
        // ------------------------------------------------------------------
        const auto clusterType = m_config.get<string>("run.vaccine_hesitancy_cluster_pooltype", "None");
        unsigned int numClustered = 0;
        unsigned int numLinked    = 0;
        vector<unsigned int> clusterIds;

        vector<bool> isLinked(pools.size(), false);
        if (clusterType != "None" && clusterType != "NA" && !clusterType.empty()) {
                const auto clusterTypeId = ContactType::ToId(clusterType);
                const auto numClusters   = m_config.get<unsigned int>("run.vaccine_hesitancy_num_clusters", 0U);
                const auto clusterRate   = m_config.get<double>("run.vaccine_hesitancy_cluster_rate", 0.0);
                if (clusterRate < 0.0 || clusterRate > 1.0) {
                        throw runtime_error("ImmunitySeeder> run.vaccine_hesitancy_cluster_rate must be in [0,1], got " +
                                            to_string(clusterRate));
                }

                if (numClusters > 0 && clusterRate > 0.0) {
                        // candidate cluster pools: non-empty pools with at least one eligible-age member
                        const auto&          clusterPools = pop->CRefPoolSys().CRefPools(clusterTypeId);
                        vector<unsigned int> candidatePools;
                        for (unsigned int i = 1; i < clusterPools.size(); i++) {
                                for (const auto& p : clusterPools[i].GetPool()) {
                                        if (isEligibleAge(p)) {
                                                candidatePools.push_back(i);
                                                break;
                                        }
                                }
                        }

                        // draw the cluster pools uniformly without replacement
                        const auto nPools = static_cast<unsigned int>(candidatePools.size());
                        const auto k      = std::min(numClusters, nPools);
                        for (unsigned int i = 0; i < k; i++) {
                                unsigned int j = i + static_cast<unsigned int>(rng.SampleUniform01() * (nPools - i));
                                if (j >= nPools) { j = nPools - 1; }
                                std::swap(candidatePools[i], candidatePools[j]);
                                clusterIds.push_back(candidatePools[i]);
                        }

                        // eligible households linked to a cluster pool, in household order
                        for (const auto poolIdx : clusterIds) {
                                for (const auto& p : clusterPools[poolIdx].GetPool()) {
                                        const auto hh = p->GetPoolId(Id::Household);
                                        if (hh < pools.size() && isEligible[hh]) {
                                                isLinked[hh] = true;
                                        }
                                }
                        }
                        vector<unsigned int> linked;
                        for (unsigned int i = 0; i < pools.size(); i++) {
                                if (isLinked[i]) {
                                        linked.push_back(i);
                                }
                        }
                        numLinked = static_cast<unsigned int>(linked.size());

                        const auto target = std::min(
                            numHesitant,
                            static_cast<unsigned int>(floor(clusterRate * static_cast<double>(numLinked) + 0.5)));
                        numClustered = markHesitant(linked, target);
                }
        }

        // the remainder: at random over all other (non-empty) households, with or without
        // children, outside the clusters
        vector<unsigned int> others;
        for (unsigned int i = 0; i < pools.size(); i++) {
                if (pools[i].size() > 0 && !isLinked[i]) {
                        others.push_back(i);
                }
        }
        const auto numRandom = markHesitant(others, numHesitant - numClustered);

        // Split into compliant and hesitant households, keeping the original order so the
        // samplers visit households in the same order as without hesitancy. Households
        // without eligible members cannot be vaccinated in this pass, so they are left out.
        compliantPools.clear();
        hesitantPools.clear();
        hesitantClusterPools.clear();
        for (unsigned int i = 0; i < pools.size(); i++) {
                if (!isEligible[i]) {
                        continue;
                }
                if (isHesitant[i] && isLinked[i]) {
                        hesitantClusterPools.push_back(pools[i]);
                } else if (isHesitant[i]) {
                        hesitantPools.push_back(pools[i]);
                } else {
                        compliantPools.push_back(pools[i]);
                }
        }

        // short report on stdout (not in the event log, to keep its format unchanged)
        unsigned int numEligibleHesitant = 0;
        for (unsigned int i = 0; i < pools.size(); i++) {
                if (isEligible[i] && isHesitant[i]) {
                        numEligibleHesitant++;
                }
        }
        std::cout << "Vaccine hesitancy: " << numClustered + numRandom << " of " << numHouseholds
                  << " households hesitant (" << numClustered << " clustered in " << clusterIds.size() << " "
                  << clusterType << " pool(s) linked to " << numLinked << " households, " << numRandom
                  << " at random); " << numEligibleHesitant << " of the " << numEligible
                  << " households with vaccine-eligible members are hesitant." << std::endl;

        return true;
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
