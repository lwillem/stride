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

#include <boost/property_tree/ptree_fwd.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>

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
        explicit Calendar(const boost::property_tree::ptree& configPt,unsigned int num_days);

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
			return m_school_closures[age][GetDayIndex(m_date)];
		}



        /// Check if distancing measures are in place for workplaces
        bool IsWorkplaceDistancingEnforced() const
        {
        	return m_workplace_distancing[GetDayIndex(m_date)] > 0.0;
		}

        /// Get distancing factor for workplaces
		double GetWorkplaceDistancingFactor() const
		{
			 return m_workplace_distancing[GetDayIndex(m_date)];
		}

        /// Check if distancing measures are in place for communities
		bool IsCommunityDistancingEnforced() const
		{
			 return m_community_distancing[GetDayIndex(m_date)] > 0.0;
		}

		/// Get distancing factor for community contacts
		double GetCommunityDistancingFactor() const
		{
  			return m_community_distancing[GetDayIndex(m_date)];
		}

		/// Check if contact tracing is place
		bool IsContactTracingActivated() const
		{
			 return m_contact_tracing[GetDayIndex(m_date)];
		}

  		/// Check if universal testing is place
		bool IsUniversalTestingActivated() const
		{
			 return m_universal_testing[GetDayIndex(m_date)];
		}

		/// Check if household clustering is allowed
		bool IsHouseholdClusteringAllowed() const
		{
			 return m_household_clustering[GetDayIndex(m_date)];
		}

		unsigned int GetNumberOfImportedCases() const
		{
			return m_imported_cases[GetDayIndex(m_date)];
		}

private:

		unsigned short int GetDayIndex(boost::gregorian::date date) const;
		unsigned short int GetDayIndex(std::string date) const;

		bool IsDatePartOfSimulation(boost::gregorian::date date) const
		{
			return m_date_start <= date && date < m_date_end;
		}

		bool IsDatePartOfSimulation(std::string date) const
		{
			return IsDatePartOfSimulation(boost::gregorian::from_simple_string(date));
		}

		/// Check if it's a public holiday.
		bool IsPublicHoliday() const
		{
			//return (std::find(m_public_holidays.begin(), m_public_holidays.end(), m_date) != m_public_holidays.end());
			return m_public_holidays[GetDayIndex(m_date)];
		}

		/// Check if it's weekend.
		bool IsWeekend() const
		{
			return (GetDayOfTheWeek() == 6 || GetDayOfTheWeek() == 0);
		}

		/// Initialize the calendar (json)
        void Initialize(const boost::property_tree::ptree& configPt);


		/// Initialize the calendar (csv)
        void Initialize_csv(const boost::property_tree::ptree& configPt);

        boost::gregorian::date              m_date;                       ///< Current simulated date.
        boost::gregorian::date              m_date_start;                 ///< Start simulation.
        boost::gregorian::date              m_date_end;                   ///< End simulation.
        std::vector<bool> m_public_holidays;            ///< Vector of public holidays
        std::vector<double> m_workplace_distancing;     ///< Vector with social distancing level enforcement at work places
        std::vector<double> m_community_distancing;     ///< Vector with social distancing level enforcement in the community
        std::vector<bool> m_contact_tracing;            ///< Vector of days with case finding measures
        std::vector<bool> m_universal_testing;          ///< Vector of days with universal testing measures
        std::vector<bool> m_household_clustering;       ///< Vector of days when household clusters are allowed

        std::vector<unsigned int>m_imported_cases; ///<Vector of days when cases are imported (~daily seeding activated)

        std::vector<std::vector<double>> m_school_closures; /// Matrix for [age x time] with social distancing at school]

};

} // namespace stride
