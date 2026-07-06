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
 *  Copyright 2020, Willem L, Kuylen E, Broeckhove J, Libin P
 */

/**
 * @file
 * Header file for the Calendar class.
 */

#pragma once

#include "contact/ContactPool.h"

#include "util/Date.h"
#include "util/Ptree.h"

#include <algorithm>
#include <cstdlib>
#include <memory>
#include <vector>


namespace stride {

/**
 * Class that keeps track of the 'state' of simulated world.
 * E.g. what day it is, holidays, quarantines, ...
 */
class Calendar
{
public:
        /// Constructor
        explicit Calendar(const stride::util::ptree& configPt,unsigned int num_days);

        /// Advance the simulated calendar by one day.
        void AdvanceDay();

        /// Current day of the month in the simulated calendar.
        std::size_t GetDay() const;

        /// Current day of the week (0 (Sunday), ..., 6 (Saturday)) in the simulated calendar.
        std::size_t GetDayOfTheWeek() const;

        /// Current month in the simulated calendar.
        std::size_t GetMonth() const;

        /// Current simulated day since the start of the simulation.
        unsigned short int GetSimulationDay() const;

        /// Current year in the simulated calendar.
        std::size_t GetYear() const;


        /// Check if today is a regular weekday (= NO weekend or holiday).
		bool IsRegularWeekday() const
		{
			return !(IsWeekend() || IsPublicHoliday());
		}

		bool IsSchoolClosed(unsigned int age) const
		{
			// if an age is not present in school closure matrix, the age is
			// not part of school ages so return false
			if(age >= m_school_closures.size()){ return false; }

			// else check calendar
			return GetSchoolDistancingFactor(age) == 1.0;
		}

		double GetSchoolDistancingFactor(unsigned int age) const
		{
			// if the age is not part of school ages, or the current day is
			// not a regular weekday, return 1 ( == full distancing)
			if(age >= m_school_closures.size() || !IsRegularWeekday()){
				return 1.0;
			}

			// else, return numerical value
			return m_school_closures[age][m_day_index];
		}


        /// Check if distancing measures are in place for workplaces
        bool IsWorkplaceDistancingEnforced() const
        {
        	return m_workplace_distancing[m_day_index] > 0.0;
		}

        /// Get distancing factor for workplaces
		double GetWorkplaceDistancingFactor() const
		{
			 return m_workplace_distancing[m_day_index];
		}

        /// Check if distancing measures are in place for communities
		bool IsCommunityDistancingEnforced() const
		{
			 return m_community_distancing[m_day_index] > 0.0;
		}

		/// Get distancing factor for community contacts
		double GetCommunityDistancingFactor() const
		{
  			return m_community_distancing[m_day_index];
		}

		/// Get distancing factor for collectivities
		double GetCollectivityDistancingFactor() const
		{
			return m_collectivity_distancing[m_day_index];
		}

		/// Get ventilation factor
		double GetVentilationFactor() const
		{
			return m_ventilation[m_day_index];
		}

		/// Check if contact tracing is in place
		bool IsContactTracingActive() const
		{
			 return m_contact_tracing[m_day_index];
		}

		/// Check if social contact survey is ongoing
		bool IsContactSurveyActive() const
		{
			 return m_contact_survey[m_day_index];
		}

		/// Check if household clustering is allowed
		bool IsHouseholdClusteringAllowed() const
		{
			 return m_household_clustering[m_day_index] > 0.0;
		}

		// Get social interaction level for household clusters
		double GetHouseholdClusteringLevel() const
		{
			 return m_household_clustering[m_day_index];
		}

		unsigned int GetNumberOfImportedCases() const
		{
			return m_imported_cases[m_day_index];
		}

		/// Update the contact reduction vectors
		void UpdateCntReduction(std::vector<double> workplace_distancing, std::vector<double> community_distancing,
                                std::vector<double> collectivity_distancing);

		void RegisterInfectedSeeds(unsigned int num_infected_seeds);

		double GetDistancingFactor(const ContactPool& pool) const;

private:

		unsigned short int GetDayIndex(util::Date date) const;
		unsigned short int GetDayIndex(std::string date) const;

		bool IsDatePartOfSimulation(util::Date date) const
		{
			return m_date_start <= date && date < m_date_end;
		}

		bool IsDatePartOfSimulation(std::string date) const
		{
			return IsDatePartOfSimulation(util::Date::FromString(date));
		}

		/// Check if it's a public holiday.
		bool IsPublicHoliday() const
		{
			//return (std::find(m_public_holidays.begin(), m_public_holidays.end(), m_date) != m_public_holidays.end());
			return m_public_holidays[m_day_index];
		}

		/// Check if it's weekend.
		bool IsWeekend() const
		{
			return (GetDayOfTheWeek() == 6 || GetDayOfTheWeek() == 0);
		}


		/// Initialize the calendar (csv)
        void Initialize_csv(const stride::util::ptree& configPt);

        util::Date              m_date;                       ///< Current simulated date.
        util::Date              m_date_start;                 ///< Start simulation.
        util::Date              m_date_end;                   ///< End simulation.
        std::vector<bool>   m_public_holidays;          ///< Vector of public holidays
        std::vector<double> m_workplace_distancing;     ///< Vector with daily social distancing level enforcement at workplaces
        std::vector<double> m_community_distancing;     ///< Vector with daily social distancing level enforcement in the community
        std::vector<double> m_collectivity_distancing;  ///< Vector with daily social distancing level enforcement in collectivities
        std::vector<bool>   m_contact_tracing;          ///< Vector of days with case finding measures
        std::vector<bool>   m_contact_survey;           ///< Vector of days to conduct a social contact survey
        std::vector<double> m_household_clustering;     ///< Vector with daily social interaction level within household clusters
		std::vector<double> m_ventilation;              ///< Vector with ventilation increase or decrease      

        std::vector<unsigned int>        m_imported_cases;  ///<Vector imported cases per day (for initial and/or daily seeding)
        std::vector<std::vector<double>> m_school_closures; /// Matrix for [age x time] with social distancing at school]

        std::size_t m_weekday;
        unsigned short int m_day;
        unsigned short int m_day_index;

};

} // namespace stride
