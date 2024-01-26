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
 * Header file for the SurveyManager class.
 */

#pragma once

#include <boost/property_tree/ptree_fwd.hpp>
#include <memory>

namespace stride {

class Population;
class Person;

namespace util {
class RnMan;
}

/**
 * Seeds the population with survey participants.
 */
class SurveyManager
{
public:
        /// Initialize Seeder.
        /// \param config         Configuration parameters.
        /// \param rnMan         Random number manager.
		SurveyManager(std::shared_ptr<Population> pop, const boost::property_tree::ptree& config, std::shared_ptr<util::RnMan> rnMan);

        /// Seeds the population with survey participants.
        /// \param pop               Population.
        void ManagePanel(std::shared_ptr<Population> pop, unsigned int simDay = 0U);

        /// Register a selected person as a survey participant
        /// \param p 				Person to register
        void RegisterParticipant(std::shared_ptr<Population> pop, Person& p, unsigned int simDay, std::string& survey_type);

        /// Register the health state of the survey participants
        void LogHealthStates();
private:
        std::shared_ptr<Population>           m_population;
		const boost::property_tree::ptree&    m_config; ///< Run config.
        std::shared_ptr<util::RnMan>          m_rn_man; ///< Random number manager.
        bool m_panel_ready;
};

} // namespace stride
