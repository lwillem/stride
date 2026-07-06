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
 *  Copyright 2018, Kuylen E, Willem L, Broeckhove J
 */

/**
 * @file
 * Implementation for the DiseaseSeeder class.
 */

#include "DiseaseSeeder.h"

#include "contact/EventLogMode.h"
#include "pop/Population.h"
#include "pop/SurveyManager.h"
#include "util/FileSys.h"
#include "util/LogUtils.h"
#include "util/StringUtils.h"


#include "util/Ptree.h"

namespace stride {

using namespace stride::ContactType;
using namespace stride::util;
using namespace std;

DiseaseSeeder::DiseaseSeeder(const ptree& config, std::shared_ptr<util::RnMan> rnMan) : m_config(config), m_rn_man(rnMan) {}

void DiseaseSeeder::ImportInfectedCases(std::shared_ptr<Population> pop, unsigned int numInfected, unsigned int simDay, const TransmissionProfile& transProfile, util::Rn& rn)
{

        // --------------------------------------------------------------
        // Add infected persons.
        // --------------------------------------------------------------
        const auto   sAgeMin     = m_config.get<double>("run.seeding_age_min", 1);
        const auto   sAgeMax     = m_config.get<double>("run.seeding_age_max", 99);
        const auto   popSize     = pop->size();
        const auto   maxPopIndex = static_cast<int>(popSize - 1);
        auto         generator   = m_rn_man->at(0U).GetUniformIntGenerator(0, maxPopIndex);
        auto&        logger      = pop->RefEventLogger();
        const EventLogMode::Id log_level   = EventLogMode::ToMode(m_config.get<string>("run.event_log_level", "None"));

        while (numInfected > 0) {
                Person& p = pop->at(static_cast<size_t>(generator()));
                if (p.GetHealth().IsSusceptible() && (p.GetAge() >= sAgeMin) && (p.GetAge() <= sAgeMax)) {
                        double rel_inf = transProfile.GetIndividualInfectiousness(rn);
                        p.GetHealth().StartInfection(p.GetId(),0,rel_inf); // TODO why is infector_id 0 here? it is logged as -1 for index cases
                        numInfected--;

                        short int startHospitalisation = -1;
						short int endHospitalisation = -1;
						if (p.GetHealth().GetStartHospitalisation()) {
							startHospitalisation = p.GetHealth().GetStartHospitalisation().value();
							endHospitalisation = p.GetHealth().GetEndHospitalisation().value();
						}

                        //TODO: make use of Infector template functions
                        if (log_level >= EventLogMode::Id::Transmissions) {
                                logger->info("[PRIM] {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}",
                                		p.GetId(), -1, p.GetAge(), -1, -1, simDay, p.GetId(),
										p.GetHealth().GetStartInfectiousness(),p.GetHealth().GetEndInfectiousness(),
										p.GetHealth().GetStartSymptomatic(),p.GetHealth().GetEndSymptomatic(),
                                        startHospitalisation,
										endHospitalisation,
                                        -1,
										p.GetHealth().GetRelativeInfectiousness(),
										p.GetHealth().GetRelativeSusceptibility(),
										0,false);

                        } else if(log_level == EventLogMode::Id::Incidence){
                        	{
								logger->info("[TRAN_M] {} {} {} {} {} {} {}",
											 p.GetAge(),
											 simDay,
											 p.GetHealth().GetStartInfectiousness(),
											 p.GetHealth().GetStartSymptomatic(),
											 p.GetHealth().GetEndSymptomatic(),
                                             startHospitalisation,
											 endHospitalisation
											 );
                        	}

                        }

                        // register as survey participant
                        //TODO: add link with logLevel
                        SurveyManager sManager(pop,m_config,m_rn_man);
                        std::string survey_type = "infection";
                        sManager.RegisterParticipant(p,simDay,survey_type);
                }
        }
}

} // namespace stride
