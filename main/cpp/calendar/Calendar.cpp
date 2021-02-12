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

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

namespace stride {

using namespace std;
using namespace boost::property_tree::json_parser;
using namespace stride::util;
using boost::property_tree::ptree;


Calendar::Calendar(const ptree& configPt,unsigned int num_days) :
		m_date(), m_date_start(), m_date_end(), m_public_holidays(num_days),
		m_workplace_distancing(num_days), m_community_distancing(num_days), m_contact_tracing(num_days),
		m_universal_testing(num_days), m_household_clustering(num_days), m_imported_cases(num_days,0U),
		m_school_closures(100, vector<double>(num_days))
{
        // Set start date
        m_date = boost::gregorian::from_simple_string(configPt.get<string>("run.start_date", "2020-01-01"));
        m_date_start = m_date;
        m_date_end = m_date + boost::gregorian::days(num_days);

        string holiday_file = configPt.get<string>("run.holidays_file", "holidays_belgium_2019_2021.csv");
        string csv_extension = "csv";

		// temporary switch
		if (IsSubstring(holiday_file, csv_extension)){
			Initialize_csv(configPt);   // csv file
		} else{
			Initialize(configPt);       // json file => default
		}
}

void Calendar::AdvanceDay()
{
        m_date = m_date + boost::gregorian::date_duration(1);
}

size_t Calendar::GetDay() const { return m_date.day(); }

size_t Calendar::GetDayOfTheWeek() const { return m_date.day_of_week(); }

size_t Calendar::GetMonth() const { return m_date.month(); }

unsigned short int Calendar::GetSimulationDay() const {

	return GetDayIndex(m_date);
}


unsigned short int Calendar::GetDayIndex(boost::gregorian::date date) const{

	if(date == m_date_start){ return 0; }

	return (date - m_date_start).days();
}

unsigned short int Calendar::GetDayIndex(std::string date) const{

	return GetDayIndex(boost::gregorian::from_simple_string(date));
}


size_t Calendar::GetYear() const { return m_date.year(); }




void Calendar::Initialize(const ptree& configPt)
{
		// include warning
		std::cout << "WARNING: JSON CALENDAR FILES WILL NOT BE SUPPORTED IN FUTURE VERSIONS... PLEASE SWITCH TO CSV" << std::endl;

		// Load json file
        ptree holidaysPt;
		const string        fName{configPt.get<string>("run.holidays_file", "holidays_flanders_2020.json")};
		const filesys::path fPath{FileSys::GetDataDir() /= fName};
		if (!is_regular_file(fPath)) {
				throw runtime_error(string(__func__) + "Holidays file " + fPath.string() + " not present.");
		}
		read_json(fPath.string(), holidaysPt);

        // Read in holidays
        for (int i = 1; i < 13; i++) {
                const auto month = to_string(i);
                const auto year  = holidaysPt.get<string>("year", "2020");
                const auto lead  = string(year).append("-").append(month).append("-");

                // read in general holidays
                for (const auto& date : holidaysPt.get_child("general." + month)) {
                        const auto d_date = string(lead).append(date.second.get_value<string>());
                         if(IsDatePartOfSimulation(d_date)){
                        	m_public_holidays[GetDayIndex(d_date)] = true;
                        }
                }

                // read in pre-school holidays
                for (const auto& date : holidaysPt.get_child("preschool." + month)) {
                        const auto d_date = string(lead).append(date.second.get_value<string>());
						if(IsDatePartOfSimulation(d_date)){
							for(int i = 0; i<= 6;i ++){
								m_school_closures[i][GetDayIndex(d_date)] = true;
							}
						}
                }
                // read in primary school holidays
				for (const auto& date : holidaysPt.get_child("primary_school." + month)) {
						const auto d_date = string(lead).append(date.second.get_value<string>());
						if(IsDatePartOfSimulation(d_date)){
							for(int i = 6; i<= 11;i ++){
								m_school_closures[i][GetDayIndex(d_date)] = true;
							}						}
				}

                // read in secondary school holidays
				for (const auto& date : holidaysPt.get_child("secondary_school." + month)) {
						const auto d_date = string(lead).append(date.second.get_value<string>());
						if(IsDatePartOfSimulation(d_date)){
							for(int i = 12; i<= 17;i ++){
								m_school_closures[i][GetDayIndex(d_date)] = true;
							}
						}
				}

				// read in college holidays
                for (const auto& date : holidaysPt.get_child("college." + month)) {
						const auto d_date = string(lead).append(date.second.get_value<string>());
						if(IsDatePartOfSimulation(d_date)){
							for(int i = 18; i<= 25;i ++){
								m_school_closures[i][GetDayIndex(d_date)] = true;
							}

						}

                }
                // read in work place distancing data (if present)
                if(holidaysPt.count("workplace_distancing") != 0){
					for (const auto& date : holidaysPt.get_child("workplace_distancing." + month)) {
							const auto d_date = string(lead).append(date.second.get_value<string>());
							if(IsDatePartOfSimulation(d_date)){
								m_workplace_distancing[GetDayIndex(d_date)] = true;
							}
					}
                }

                // read in community distancing data (if present)
				if(holidaysPt.count("community_distancing") != 0){
					for (const auto& date : holidaysPt.get_child("community_distancing." + month)) {
							const auto d_date = string(lead).append(date.second.get_value<string>());
							if(IsDatePartOfSimulation(d_date)){
								m_community_distancing[GetDayIndex(d_date)] = true;
							}

					}
				}

				// read in contact tracing data (if present)
				if(holidaysPt.count("contact_tracing") != 0){
					for (const auto& date : holidaysPt.get_child("contact_tracing." + month)) {
							const auto d_date = string(lead).append(date.second.get_value<string>());
							if(IsDatePartOfSimulation(d_date)){
								m_contact_tracing[GetDayIndex(d_date)] = true;
							}

					}
				}

				// read in universal testing data (if present)
				if(holidaysPt.count("universal_testing") != 0){
					for (const auto& date : holidaysPt.get_child("universal_testing." + month)) {
							const auto d_date = string(lead).append(date.second.get_value<string>());
							if(IsDatePartOfSimulation(d_date)){
								m_universal_testing[GetDayIndex(d_date)] = true;
							}

					}
				}

				// read in household clustering data (if present)
				if(holidaysPt.count("household_clustering") != 0){
					for (const auto& date : holidaysPt.get_child("household_clustering." + month)) {
							const auto d_date = string(lead).append(date.second.get_value<string>());
							if(IsDatePartOfSimulation(d_date)){
								m_household_clustering[GetDayIndex(d_date)] = true;
							}

					}
				}

				// read in imported cases
				if(holidaysPt.count("import_cases") != 0){
					for (const auto& date : holidaysPt.get_child("import_cases." + month)) {
						const auto d_date = string(lead).append(date.second.get_value<string>());
						unsigned int num_cases = configPt.get<unsigned int>("run.num_daily_imported_cases",0);
						if(IsDatePartOfSimulation(d_date)){
							m_imported_cases[GetDayIndex(d_date)] = num_cases;
						}
					}
				} else { // if no calendar info present, use the same value throughout the simulation
					unsigned int num_cases = configPt.get<unsigned int>("run.num_daily_imported_cases",0);
					for (unsigned int day_index = 0 ; day_index < m_imported_cases.size() ; day_index++){
						m_imported_cases[day_index] = num_cases;
					}
				}
        }
}

void Calendar::Initialize_csv(const ptree& configPt)
{
        // Load csv file
		const auto fileName = configPt.get<string>("run.holidays_file", "holidays_belgium_2019_2021.csv");
		const filesys::path filePath{FileSys::GetDataDir() /= fileName};
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
				const auto date = boost::gregorian::from_simple_string(date_str);

				// check date
				if(IsDatePartOfSimulation(date_str)){

					// convert value into boolean
					const bool value_boolean = value == 1.0;

					if(category == "general")              {  m_public_holidays[GetDayIndex(date)] = value_boolean; }
					if(category == "schools_closed")       {  m_school_closures[age][GetDayIndex(date)] = value; }
					if(category == "workplace_distancing") {  m_workplace_distancing[GetDayIndex(date)] = value; }
					if(category == "community_distancing") {  m_community_distancing[GetDayIndex(date)] = value; }
					if(category == "household_clustering") {  m_household_clustering[GetDayIndex(date)] = value_boolean; }
					if(category == "contact_tracing")      {  m_contact_tracing[GetDayIndex(date)] = value_boolean; }
					if(category == "universal_testing")    {  m_universal_testing[GetDayIndex(date)] = value_boolean; }
					if(category == "imported_cases")
					{
						unsigned int num_cases = configPt.get<unsigned int>("run.num_daily_imported_cases",0);
						m_imported_cases[GetDayIndex(date)] = num_cases;
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


} // namespace stride
