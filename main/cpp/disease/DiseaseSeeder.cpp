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
#include "pop/SurveySeeder.h"
#include "util/FileSys.h"
#include "util/LogUtils.h"
#include "util/StringUtils.h"


#include <boost/property_tree/ptree.hpp>

namespace stride {

using namespace boost::property_tree;
using namespace stride::ContactType;
using namespace stride::util;
using namespace std;

DiseaseSeeder::DiseaseSeeder(const ptree& config, RnMan& rnMan) : m_config(config), m_rn_man(rnMan) {}

void DiseaseSeeder::Seed(std::shared_ptr<Population> pop, const TransmissionProfile& transProfile, util::RnHandler& rnHandler)
{

        // --------------------------------------------------------------
        // Get number of infected persons to seed
        // --------------------------------------------------------------
        const auto   sRate       = m_config.get<double>("run.seeding_rate",0);
        const auto   popSize     = pop->size();
        auto numInfected = static_cast<unsigned int>(floor(static_cast<double>(popSize) * sRate));


        // option to provide a number of cases, instead of a seeding rate.
	    const auto numInfectedSeeds = m_config.get<unsigned int>("run.num_infected_seeds",0);
	    if(numInfectedSeeds > 0){
	    	numInfected = numInfectedSeeds;
	    }

	    //------------------------------------------------
		// Check compatibility options issues
		//------------------------------------------------
	    if(sRate > 0){
			std::cout << "The 'seeding_rate' parameter is deprecated. Please use 'num_infected_seeds'" << std::endl;
			std::cerr << "The 'seeding_rate' parameter is deprecated. Please use 'num_infected_seeds'" << std::endl;
		}
	    if(numInfectedSeeds > 0 && sRate > 0){
    		std::cout << "The 'seeding_rate' parameter is overruled by 'num_infected_seeds'" << std::endl;
    		std::cerr << "The 'seeding_rate' parameter is overruled by 'num_infected_seeds' " << std::endl;
    	}

	    //------------------------------------------------
		// Check validity of input data.
		//------------------------------------------------
	    if (numInfected > popSize) {
				throw runtime_error(string(__func__) + "> Bad input data for infected seeding an/or initially infected cases.");
		}

        // --------------------------------------------------------------
        // Add infected seeds to the population
        // --------------------------------------------------------------
        ImportInfectedCases(pop, numInfected, 0, transProfile, rnHandler);
}

void DiseaseSeeder::ImportInfectedCases(std::shared_ptr<Population> pop, unsigned int numInfected, unsigned int simDay, const TransmissionProfile& transProfile, util::RnHandler& rnHandler)
{

        // --------------------------------------------------------------
        // Add infected persons.
        // --------------------------------------------------------------
        const auto   sAgeMin     = m_config.get<double>("run.seeding_age_min", 1);
        const auto   sAgeMax     = m_config.get<double>("run.seeding_age_max", 99);
        const auto   popSize     = pop->size();
        const auto   maxPopIndex = static_cast<int>(popSize - 1);
        auto         generator   = m_rn_man.GetUniformIntGenerator(0, maxPopIndex, 0U);
        auto&        logger      = pop->RefEventLogger();
        const EventLogMode::Id log_level   = EventLogMode::ToMode(m_config.get<string>("run.event_log_level", "None"));

        while (numInfected > 0) {
                Person& p = pop->at(static_cast<size_t>(generator()));
                if (p.GetHealth().IsSusceptible() && (p.GetAge() >= sAgeMin) && (p.GetAge() <= sAgeMax)) {
                        double rel_inf = transProfile.GetIndividualInfectiousness(rnHandler);
                        p.GetHealth().StartInfection(p.GetId(),0,rel_inf);
                        numInfected--;

                        //TODO: make use of Infector template functions
                        if (log_level >= EventLogMode::Id::Transmissions) {
                                logger->info("[PRIM] {} {} {} {} {} {} {} {} {} {} {} {} {} {}",
                                		p.GetId(), -1, p.GetAge(), -1, -1, simDay, p.GetId(),
										p.GetHealth().GetStartInfectiousness(),p.GetHealth().GetEndInfectiousness(),
										p.GetHealth().GetStartSymptomatic(),p.GetHealth().GetEndSymptomatic(), -1,
										p.GetHealth().GetRelativeInfectiousness(),
										p.GetHealth().GetRelativeSusceptibility());

                        } else if(log_level == EventLogMode::Id::Incidence){
                        	{
								logger->info("[TRAN_M] {} {} {} {} {}",
											 p.GetAge(),
											 simDay,
											 p.GetHealth().GetStartInfectiousness(),
											 p.GetHealth().GetStartSymptomatic(),
											 p.GetHealth().GetEndSymptomatic()
											 );
							}
                        }

                        // register as survey participant
                        //TODO: add link with logLevel
                        SurveySeeder sSeeder(m_config,m_rn_man);
                        sSeeder.RegisterParticipant(pop,p);
                }
        }
}


} // namespace stride
