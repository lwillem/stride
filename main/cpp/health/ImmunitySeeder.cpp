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

        // Mark vaccine-hesitant households before the vaccine pass. Natural immunity
        // above is unaffected -- hesitancy only governs whether a household's children
        // are eligible for vaccination.
        MarkHesitantHouseholds(pop);

        const auto vaccinationProfile = m_config.get<std::string>("run.vaccine_profile");
        if(vaccinationProfile == "Teachers"){
        	Vaccinate("vaccine", "Random", pop->CRefPoolSys().CRefPools<Id::School>(),pop);
        } else {
        	Vaccinate("vaccine", vaccinationProfile, pop->CRefPoolSys().CRefPools<Id::Household>(),pop);
        }
}

void ImmunitySeeder::MarkHesitantHouseholds(std::shared_ptr<Population> pop)
{
        // Optional config key; if absent or zero, no households are marked and
        // behaviour is identical to before this change.
        const auto hesitancyRate = m_config.get<double>("run.vaccine_hesitancy_rate", 0.0);
        if (hesitancyRate <= 0.0) {
                return;
        }

        for (auto& hh : pop->CRefPoolSys().CRefPools<Id::Household>()) {
                if (m_rn_man->at(0U).SampleUniform01() < hesitancyRate) {
                        m_hesitant_households.insert(hh.GetId());
                }
        }
}

void ImmunitySeeder::Vaccinate(const std::string& immunityType, const std::string& immunizationProfile,
                              const SegmentedVector<ContactPool>& immunityPools,std::shared_ptr<Population> pop)
{
        std::vector<double> immunityDistribution;
        double              linkProbability = 0;

        // retrieve the maximum age in the population
        unsigned int maxAge = pop->GetMaxAge();

        // Exclude vaccine-hesitant households from the vaccine pass only. Household-type
        // pools are the only ones checked against m_hesitant_households (the "Teachers"
        // path passes School pools, which are left untouched).
        SegmentedVector<ContactPool> eligibleImmunityPools;
        for (auto& c : immunityPools) {
                if (immunityType == "vaccine" && c.GetType() == Id::Household &&
                    m_hesitant_households.count(c.GetId()) > 0) {
                        continue;
                }
                eligibleImmunityPools.push_back(c);
        }

        if (immunizationProfile == "AgeDependent") {
                        const auto   immunityFile = m_config.get<string>("run." + ToLower(immunityType) + "_distribution_file");
                        const ptree& immunity_pt  = FileSys::ReadPtreeFile(immunityFile);

                        linkProbability = m_config.get<double>("run." + ToLower(immunityType) + "_link_probability");

                        for (unsigned int index_age = 0; index_age <= maxAge; index_age++) {
                                auto immunityRate = immunity_pt.get<double>("immunity.age" + std::to_string(index_age));
                                immunityDistribution.push_back(immunityRate);
                        }
                        Random(eligibleImmunityPools, immunityDistribution, linkProbability, pop, false); // CHANGED: eligibleImmunityPools

		} else if(immunizationProfile == "Random" || immunizationProfile == "Cocoon") {

			// Initialize new ContactPool vector
			SegmentedVector<ContactPool> immunityPools_selection;

			// immunizationProfile == Random: copy all contact pools
			// immunizationProfile == Cocoon: copy all contact pools with an infant
			for (auto& c : eligibleImmunityPools) { // CHANGED: iterate eligibleImmunityPools instead of immunityPools
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

        // Previously a random pool was drawn WITH replacement on every step via
        // GetUniformIntGenerator(), which meant popular pools got redrawn and fully
        // reshuffled many times as quotas neared their target (coupon-collector cost --
        // this is the main cause of slowness at 90%+ rates). Instead, visit households
        // in "laps": each lap is one shuffled pass over ALL households, so no household
        // is revisited until every other household has had a turn in that lap. Multiple
        // laps still happen (needed when immunityLinkProbability is low, since each
        // household visit may stop after just one member), but coverage per lap is now
        // uniform instead of random-with-replacement, which sharply cuts the number of
        // laps needed. Clustering behaviour (immunityLinkProbability) is unchanged.
        const auto poolsSize = static_cast<int>(pools.size());
        vector<unsigned int> poolOrder(poolsSize);
        iota(poolOrder.begin(), poolOrder.end(), 0U);

        while (numImmune > 0) {
                m_rn_man->at(0U).Shuffle(poolOrder);

                for (unsigned int idx = 0; idx < poolOrder.size() && numImmune > 0; idx++) {
                        const ContactPool&   p_pool = pools[poolOrder[idx]];
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
                                // random draw to continue in this pool or to move to the next one in this lap
                                if (m_rn_man->at(0).SampleUniform01() < (1 - immunityLinkProbability)) {
                                        break;
                                }
                        }
                }
        }
}


} // namespace stride