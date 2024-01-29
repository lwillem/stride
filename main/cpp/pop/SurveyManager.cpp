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
 * Implementation for the SurveyManager class.
 */

#include "SurveyManager.h"

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

SurveyManager::SurveyManager(std::shared_ptr<Population> pop, const ptree& config, std::shared_ptr<RnMan> rnMan) :
		m_population(pop), m_rn_man(rnMan), m_panel_ready(false){

	m_resample_panel   = config.get<unsigned int>("run.contact_survey_resample",0) == 1;
	m_is_survey_active = EventLogMode::ToMode(config.get<string>("run.event_log_level", "None")) != EventLogMode::Id::None;
	m_num_participants = config.get<unsigned int>("run.num_participants_survey");

	m_quota_symptomatic = config.get<double>("run.contact_survey_quota_symptomatic",0);
	m_quota_symptomatic = m_quota_symptomatic > 1 ? 1 : m_quota_symptomatic;
	m_quota_symptomatic = m_quota_symptomatic < 0 ? 0 : m_quota_symptomatic;

}

void SurveyManager::ManagePanel(unsigned int simDay)
{
	if(m_is_survey_active){
		if(!m_panel_ready) {
			m_panel_ready = !m_resample_panel;
			SampleParticipantsWithQuota(simDay);
		}
		LogHealthStates(simDay);
	}
}

void SurveyManager::SampleParticipantsWithQuota(unsigned int simDay){

	// Clear panel (if existing)
	ClearPanel();

	Population& population  = *m_population;

	const auto  popCount    = static_cast<unsigned int>(population.size() - 1);
	auto  numSurveyed       = m_num_participants;

	assert((popCount >= 1U) && "SurveySeeder> Population count zero unacceptable.");
	assert((popCount >= numSurveyed) && "SurveySeeder> Pop count has to exceed the number of survey participants.");

	// Make sure the number of survey participants does not outnumber the population size (else no survey)
	if(popCount < numSurveyed){
		numSurveyed = 1;
	}

	unsigned int numSymptomatic     = static_cast<unsigned int>(std::round(numSurveyed * m_quota_symptomatic));
	if(numSymptomatic > population.CountSymptomaticCases()){
			numSymptomatic = population.CountSymptomaticCases();
	}
	unsigned int numNonSymptomatic = numSurveyed - numSymptomatic;

	// loop over population in random order
	vector<unsigned int> indices(popCount);
	iota(indices.begin(), indices.end(), 0U);
	m_rn_man->at(0U).Shuffle(indices);

	for (unsigned int i_p = 0; i_p < popCount && numSurveyed > 0; i_p++) {
		Person& p = population[indices[i_p]];

		if(m_quota_symptomatic > 0 && p.GetHealth().IsSymptomatic()){
			if(numSymptomatic > 0) { numSymptomatic--;} else {continue;}
		} else {
			if(numNonSymptomatic > 0) { numNonSymptomatic--;} else {continue;}
		}

		// register new participant
		std::string survey_type = "contacts";
		RegisterParticipant(p,simDay,survey_type);

		// update number of remaining samples
		numSurveyed--;

	}
}

void SurveyManager::RegisterParticipant(Person& p, unsigned int simDay ,std::string& survey_type)
{

	if (m_is_survey_active) {

		Population& population  = *m_population;

		auto&       poolSys     = population.CRefPoolSys();
		auto&       logger      = population.RefEventLogger();

		// set person flag to be survey participant
		p.ParticipateInSurvey();

		// log person details
		const auto h    = p.GetHealth();
		const auto pHH  = p.GetPoolId(Id::Household);
		const auto pS   = p.GetPoolId(Id::School);
		const auto pW   = p.GetPoolId(Id::Workplace);
		const auto pPC  = p.GetPoolId(Id::CommunityWeekend);
		const auto pSC  = p.GetPoolId(Id::CommunityWeekday);
		const auto pHC  = p.GetPoolId(Id::HouseholdCluster);
		const auto pCol = p.GetPoolId(Id::Collectivity);

		logger->info("[PART] {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}",
			 p.GetId(), p.GetAge(), pHH, pS, pW, pHC, pCol, h.IsSusceptible(), h.IsInfected(), h.IsInfectious(),
			 h.IsSymptomatic(),
			 h.IsRecovered(), p.IsImmune(), h.GetStartInfectiousness(), h.GetStartSymptomatic(),
			 h.GetStartHospitalisationValue(),
			 h.GetEndInfectiousness(), h.GetEndSymptomatic(),
			 h.GetEndHospitalisationValue(),
			 poolSys.CRefPools<Id::Household>()[pHH].GetPool().size(),
			 poolSys.CRefPools<Id::School>()[pS].GetPool().size(),
			 poolSys.CRefPools<Id::Workplace>()[pW].GetPool().size(),
			 poolSys.CRefPools<Id::CommunityWeekend>()[pPC].GetPool().size(),
			 poolSys.CRefPools<Id::CommunityWeekday>()[pSC].GetPool().size(),
			 simDay, survey_type
			 );
	 }
}


//TODO: check omp options
void SurveyManager::ClearPanel(){

	Population& population  = *m_population;

	for (size_t i = 0; i < population.size(); ++i) {
    	population[i].QuitSurvey();
	}

}

//TODO: check omp options
void SurveyManager::LogHealthStates(unsigned int simDay){

	if (m_is_survey_active) {

		Population& population  = *m_population;
		auto&       logger      = population.RefEventLogger();

		for (auto& p : population) {
			if(p.IsSurveyParticipant()){
				// log person details
				const auto h    = p.GetHealth();
				logger->info("[HEALTH] {} {} {} {} {} {} {} {}",
								 p.GetId(), simDay, h.IsSusceptible(), h.IsInfected(), h.IsInfectious(),
								 h.IsSymptomatic(), h.IsRecovered(), p.IsImmune()
								 );
					}

			}
		}
}


} // namespace stride
