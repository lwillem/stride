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
 * Implementation file for the Calendar class.
 */

#include "Calendar.h"

#include "util/FileSys.h"
#include "util/StringUtils.h"

#include "util/Ptree.h"
#include <fstream>

namespace stride {

using namespace std;
using namespace stride::util;
using stride::util::ptree;


Calendar::Calendar(const ptree& configPt,unsigned int num_days) :
		m_date(), m_date_start(), m_date_end(), m_public_holidays(num_days),
		m_workplace_distancing(num_days), m_community_distancing(num_days), m_collectivity_distancing(num_days),
		m_contact_tracing(num_days),
		m_contact_survey(num_days), m_household_clustering(num_days), m_ventilation(num_days),
		m_imported_cases(num_days,0U),
		m_school_closures(100, vector<double>(num_days)),
        //
        m_weekday(), m_day(), m_day_index()
{
        // Set start date
        m_date = util::Date::FromString(configPt.get<string>("run.start_date", "2020-01-01"));
        m_date_start = m_date;
        m_date_end = m_date + static_cast<int>(num_days);
        //
        m_weekday = m_date.DayOfWeek();
        m_day = 0;
        m_day_index = GetDayIndex(m_date);

        string holiday_file = configPt.get<string>("run.holidays_file", "holidays_belgium_2019_2021.csv");
        Initialize_csv(configPt);   // csv file
}

void Calendar::AdvanceDay()
{
        m_date = m_date + 1;
        m_weekday = m_date.DayOfWeek();
        m_day = m_day + 1;
        m_day_index = GetDayIndex(m_date); // TODO: m_day_index = m_day, remove?
}

size_t Calendar::GetDay() const { return m_date.Day(); }

size_t Calendar::GetDayOfTheWeek() const { return m_weekday; }

size_t Calendar::GetMonth() const { return m_date.Month(); }

unsigned short int Calendar::GetSimulationDay() const {

	return m_day;
}


unsigned short int Calendar::GetDayIndex(util::Date date) const{

	return static_cast<unsigned short int>(date - m_date_start);
}

unsigned short int Calendar::GetDayIndex(std::string date) const{

	return GetDayIndex(util::Date::FromString(date));
}


size_t Calendar::GetYear() const { return m_date.Year(); }



void Calendar::Initialize_csv(const ptree& configPt)
{
        // Load csv file
		const auto fileName = configPt.get<string>("run.holidays_file", "data/holidays_belgium_2019_2021.csv");
		const filesys::path filePath{fileName};
		if (!is_regular_file(filePath)) {
				throw runtime_error(string(__func__) + "> Holidays file " + filePath.string() + " not present.");
		}

        ifstream calendarFile;
		calendarFile.open(filePath.string());
		if (!calendarFile.is_open()) {
				throw runtime_error(string(__func__) + "> Error when opening calendar file " + filePath.string());
		}

		// do we need to add "imported cases" later on?
		bool bool_no_imported_cases_dates = true;

		string line;
		getline(calendarFile, line); // step over file header

		// set and check line separator
		string line_sep = ",";
		if(!IsSubstring(line, line_sep)){
			throw runtime_error(string(__func__) + "> Error when parsing calendar file " + filePath.string() + ": no separator ',' present");
		}

		while (getline(calendarFile, line)) {

				const auto calendar_item        = Split(line, line_sep);
				const auto category             = FromString<string>(calendar_item[0]);
				const auto date_str             = FromString<string>(calendar_item[1]);
				const double value              = FromString<double>(calendar_item[2]);
				//const auto type                 = FromString<string>(calendar_item[3]);
				const auto age                  = FromString<unsigned int>(calendar_item[4]);

				// convert date
				const auto date = util::Date::FromString(date_str);
                const auto date_index = GetDayIndex(date);

				// check date
				if(IsDatePartOfSimulation(date_str)){

					// convert value into boolean
					const bool value_boolean = value == 1.0;

					if(category == "general")              {  m_public_holidays[date_index] = value_boolean; }
					if(category == "schools_closed")       {  m_school_closures[age][date_index] = value; }
					if(category == "workplace_distancing") {  m_workplace_distancing[date_index] = value; }
					if(category == "community_distancing") {  m_community_distancing[date_index] = value; }
					if(category == "collectivity_distancing"){m_collectivity_distancing[date_index] = value; }
					if(category == "household_clustering") {  m_household_clustering[date_index] = value; }
					if(category == "contact_tracing")      {  m_contact_tracing[date_index] = value_boolean; }
					if(category == "ventilation")          {  m_ventilation[date_index] = value; }
					if(category == "contact_survey")       {  m_contact_survey[date_index] = value_boolean;  }
					if(category == "imported_cases")
					{
						unsigned int num_cases = FromString<unsigned int>(calendar_item[2]);
						m_imported_cases[date_index] = num_cases;
					}

				} // end if valid date

				// check if "imported cases" is present in the calendar file
				if(category == "imported_cases"){
					bool_no_imported_cases_dates = false;
				}
			} // end iteration over all lines


		// special case if "imported cases" is not present in the calendar file
		if(bool_no_imported_cases_dates){
			unsigned int num_cases = configPt.get<unsigned int>("run.num_daily_imported_cases",0);
			for (unsigned int day_index = 0 ; day_index < m_imported_cases.size() ; day_index++){
				m_imported_cases[day_index] = num_cases;
			}
		}

		// close file stream
		calendarFile.close();
}

void Calendar::RegisterInfectedSeeds(unsigned int num_infected_seeds) {
	m_imported_cases[0] = num_infected_seeds;
}

//double Calendar::GetDistancingFactor(ContactType::Id typ) {
double Calendar::GetDistancingFactor(const ContactPool& pool) const {

	ContactType::Id cType = pool.GetType();

	double typ_distancing_factor = 0;
	if (cType == ContactType::Id::Workplace) {
		// account for physical distancing at work
		typ_distancing_factor = GetWorkplaceDistancingFactor();
	} else if (cType == ContactType::Id::CommunityWeekend ||
				cType == ContactType::Id::CommunityWeekday) {
		// account for physical distancing in the community
		typ_distancing_factor = GetCommunityDistancingFactor();
	} else if (cType == ContactType::Id::School) {
		// account for physical distancing at school
		typ_distancing_factor = GetSchoolDistancingFactor(pool.GetMinAge());
	} else if (cType == ContactType::Id::Collectivity) {
		// account for physical distancing in the collectivity
		typ_distancing_factor = GetCollectivityDistancingFactor();
	} else if (cType == ContactType::Id::HouseholdCluster) {
		// account for contact intensity in household clusters
		typ_distancing_factor = 1-GetHouseholdClusteringLevel();
	}

	return(typ_distancing_factor);
}

void Calendar::UpdateCntReduction(std::vector<double> workplace_distancing, std::vector<double> community_distancing,
                                  std::vector<double> collectivity_distancing)
{
        m_workplace_distancing = workplace_distancing;
        m_community_distancing = community_distancing;
        m_collectivity_distancing = collectivity_distancing;
}


} // namespace stride
