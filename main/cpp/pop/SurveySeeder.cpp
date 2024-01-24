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
 *  Copyright 2024
 */

/**
 * @file
 * Implementation for the SurveySeeder class.
 */

#include "SurveySeeder.h"

#include "contact/EventLogMode.h"
#include "pop/Population.h"
#include "util/Exception.h"
#include "util/RnMan.h"
#include "util/StringUtils.h"

#include <boost/property_tree/ptree.hpp>
#include <cassert>

using namespace boost::property_tree;
using namespace stride::util;
using namespace stride::ContactType;
using namespace std;

namespace stride {

SurveySeeder::SurveySeeder(const ptree& config, std::shared_ptr<RnMan> rnMan) : m_config(config), m_rn_man(rnMan) {}

shared_ptr<Population> SurveySeeder::Seed(shared_ptr<Population> pop)
{
	const EventLogMode::Id logLevel   = EventLogMode::ToMode(m_config.get<string>("run.event_log_level", "None"));
	if (logLevel != EventLogMode::Id::None) {
		Population& population  = *pop;
		const auto  popCount    = static_cast<unsigned int>(population.size() - 1);
		auto  numSurveyed = m_config.get<unsigned int>("run.num_participants_survey");

		assert((popCount >= 1U) && "SurveySeeder> Population count zero unacceptable.");
		assert((popCount >= numSurveyed) && "SurveySeeder> Pop count has to exceed the number of survey participants.");

		// Make sure the number of survey participants does not outnumber the population size (else no survey)
		if(popCount < numSurveyed){
			numSurveyed = 1;

		}

		// Use while-loop to get 'participants' unique participants (default sampling is with replacement).
		// A for loop will not do because we might draw the same person twice.
		auto numSamples = 0U;
		auto generator  = m_rn_man->GetUniformIntGenerator(0, static_cast<int>(popCount), 0U);

		while (numSamples < numSurveyed) {
				Person& p = population[generator()];
				if (p.IsSurveyParticipant()) {
						continue;
				}

				// register new participant
				std::string survey_type = "contacts";
				RegisterParticipant(pop,p,survey_type);

				// update number of remaining samples
				numSamples++;
		}
	}
	return pop;
}

void SurveySeeder::RegisterParticipant(std::shared_ptr<Population> pop, Person& p, std::string& survey_type)
{

	const EventLogMode::Id logLevel   = EventLogMode::ToMode(m_config.get<string>("run.event_log_level", "None"));
	if (logLevel != EventLogMode::Id::None) {
		Population& population  = *pop;

		auto&       poolSys     = population.CRefPoolSys();
		auto&       logger      = population.RefEventLogger();

		// set person flag to be survey participant
		p.ParticipateInSurvey();

		// log person details
		const auto h    = p.GetHealth();
		const auto pHH  = p.GetPoolId(Id::Household);
		const auto pS   = p.GetPoolId(Id::School);
		const auto pW   = p.GetPoolId(Id::Workplace);
		const auto pPC  = p.GetPoolId(Id::PrimaryCommunity);
		const auto pSC  = p.GetPoolId(Id::SecondaryCommunity);
		const auto pHC  = p.GetPoolId(Id::HouseholdCluster);
		const auto pCol = p.GetPoolId(Id::Collectivity);

		logger->info("[PART] {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}",
			 p.GetId(), p.GetAge(), pHH, pS, pW, pHC, pCol, h.IsSusceptible(), h.IsInfected(), h.IsInfectious(),
			 h.IsRecovered(), p.IsImmune(), h.GetStartInfectiousness(), h.GetStartSymptomatic(),
			 h.GetStartHospitalisationValue(),
			 h.GetEndInfectiousness(), h.GetEndSymptomatic(),
			 h.GetEndHospitalisationValue(),
			 poolSys.CRefPools<Id::Household>()[pHH].GetPool().size(),
			 poolSys.CRefPools<Id::School>()[pS].GetPool().size(),
			 poolSys.CRefPools<Id::Workplace>()[pW].GetPool().size(),
			 poolSys.CRefPools<Id::PrimaryCommunity>()[pPC].GetPool().size(),
			 poolSys.CRefPools<Id::SecondaryCommunity>()[pSC].GetPool().size(),
			 survey_type
			 );
	 }
}


} // namespace stride
