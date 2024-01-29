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
		/// Initialize SurveyManager.
        /// \param config         Configuration parameters.
        /// \param rnMan         Random number manager.
		SurveyManager(std::shared_ptr<Population> pop, const boost::property_tree::ptree& config, std::shared_ptr<util::RnMan> rnMan);

        /// Manage the survey participants in the given population.
        /// \param pop               Population.
        void ManagePanel(unsigned int simDay = 0U);

        /// Clear survey participant panel
        /// \param pop 				Population
        void ClearPanel();

        void SampleParticipants(unsigned int simDay);

        /// Register a selected person as a survey participant
		/// \param p 				Person to register
		void RegisterParticipant(Person& p, unsigned int simDay, std::string& survey_type);

		/// Register the health state of the survey participants
        void LogHealthStates(unsigned int simDay);


private:
        std::shared_ptr<Population>           m_population;       ///< Link to the population object
        std::shared_ptr<util::RnMan>          m_rn_man;           ///< Random number manager.
        bool                                  m_panel_ready;      ///< Is the panel ready?
        bool                                  m_resample_panel;   ///< Is a new panel for each survey round needed?
        bool                                  m_is_survey_active; ///< Does the log level allow survey activities
        unsigned int                          m_num_participants; ///< Number of participants in the survey
};

} // namespace stride
