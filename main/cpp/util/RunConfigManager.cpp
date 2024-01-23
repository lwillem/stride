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
 *  Copyright 2020, Kuylen E, Willem L, Broeckhove J
 */
/**
 * @file
 * Implementation of ptree utils.
 */

#include "RunConfigManager.h"

#include "util/ConfigInfo.h"
#include "util/FileSys.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <initializer_list>
#include <map>
#include <sstream>
#include <string>

using namespace boost::property_tree;
using namespace boost::property_tree::xml_parser;
using namespace std;

namespace stride {
namespace util {

ptree RunConfigManager::Create(const std::string& configName)
{
        static map<string, string (*)()> creators{
            // clang-format off
                make_pair("TestsInfluenza",  &CreateTestsInfluenza),
                make_pair("TestsMeasles",    &CreateTestsMeasles),
				make_pair("TestsCovid19",    &CreateTestsCovid19)
        }; // clang-format on
        return FromString(creators.at(configName)());
}


vector<unsigned int> RunConfigManager::CreateNumThreads(unsigned int maxNum)
{
        maxNum = max(maxNum, ConfigInfo::NumberAvailableThreads());
        vector<unsigned int> num{1U};

        if (4 > maxNum && maxNum >= 2) {
                num.push_back(2U);
        }
        if (8 > maxNum && maxNum >= 4) {
                num.push_back(2U);
                num.push_back(4U);
        }
        if (12 > maxNum && maxNum >= 8) {
                num.push_back(2U);
                num.push_back(4U);
                num.push_back(8U);
        }
        if (16 > maxNum && maxNum >= 12) {
                num.push_back(4U);
                num.push_back(8U);
                num.push_back(12U);
        }
        if (maxNum >= 16) {
                num.push_back(4U);
                num.push_back(8U);
                num.push_back(12U);
                num.push_back(16U);
        }

        return num;
}

string RunConfigManager::CreateTestsInfluenza()
{
        return R"###(
<?xml version="1.0" encoding="utf-8"?>
<run>
        <age_contact_matrix_file>contact_matrix_flanders_conditional_teachers.xml</age_contact_matrix_file>
        <event_log_level>None</event_log_level>
        <event_output_file>false</event_output_file>
        <disease_config_file>disease_influenza.xml</disease_config_file>
        <holidays_file>holidays_none.csv</holidays_file>
        <immunity_link_probability>0</immunity_link_probability>
        <immunity_profile>Random</immunity_profile>
        <immunity_rate>0</immunity_rate>
        <num_days>30</num_days>
        <num_participants_survey>10</num_participants_survey>
        <num_threads>1</num_threads>
        <output_prefix></output_prefix>
        <population_file>pop_belgium600k_c500_teachers_censushh.csv</population_file>
        <population_type>default</population_type>
        <rng_seed>4730214</rng_seed>
        <r0>3</r0>
        <num_infected_seeds>540</num_infected_seeds>
        <seeding_age_min>1</seeding_age_min>
        <seeding_age_max>99</seeding_age_max>
        <start_date>2020-04-01</start_date>
        <stride_log_level>info</stride_log_level>
        <track_index_case>false</track_index_case>
        <use_install_dirs>true</use_install_dirs>
        <vaccine_profile>None</vaccine_profile>
</run>
        )###";
}

string RunConfigManager::CreateTestsMeasles()
{
        return R"###(
<?xml version="1.0" encoding="utf-8"?>
<run>
        <age_contact_matrix_file>contact_matrix_flanders_conditional_teachers.xml</age_contact_matrix_file>
        <event_log_level>Transmissions</event_log_level>
        <event_output_file>false</event_output_file>
        <disease_config_file>disease_measles_adaptive_behavior.xml</disease_config_file>
        <holidays_file>holidays_belgium_2019_2021.csv</holidays_file>
        <immunity_rate>0.80000000000000004</immunity_rate>
        <immunity_profile>None</immunity_profile>
        <num_days>50</num_days>
        <num_participants_survey>10</num_participants_survey>
        <num_threads>1</num_threads>
        <output_prefix></output_prefix>
        <population_file>pop_belgium600k_c500_teachers_censushh.csv</population_file>
        <population_type>default</population_type>
        <rng_seed>1097253,2387652,9963540,4730214</rng_seed>
        <seeding_age_max>99</seeding_age_max>
        <seeding_age_min>1</seeding_age_min>
        <num_infected_seeds>1200</num_infected_seeds>
        <start_date>2020-03-29</start_date>
        <stride_log_level>info</stride_log_level>
        <track_index_case>false</track_index_case>
        <use_install_dirs>true</use_install_dirs>
        <vaccine_link_probability>0</vaccine_link_probability>
        <vaccine_profile>Random</vaccine_profile>
        <vaccine_rate>0.80000000000000004</vaccine_rate>
</run>
        )###";
}

string RunConfigManager::CreateTestsCovid19()
{
        return R"###(
<?xml version="1.0" encoding="utf-8"?>
<run>
        <age_contact_matrix_file>contact_matrix_flanders_conditional_teachers.xml</age_contact_matrix_file>
        <event_log_level>None</event_log_level>
        <event_output_file>false</event_output_file>
        <disease_config_file>disease_covid19_age.xml</disease_config_file>
        <holidays_file>holidays_belgium_2019_2021.csv</holidays_file>
        <immunity_link_probability>0</immunity_link_probability>
        <immunity_profile>Random</immunity_profile>
        <immunity_rate>0</immunity_rate>
        <num_days>31</num_days>
        <num_participants_survey>10</num_participants_survey>
        <num_threads>1</num_threads>
        <output_prefix></output_prefix>
        <population_file>pop_belgium600k_c500_teachers_censushh.csv</population_file>
        <population_type>default</population_type>
        <rng_seed>4730214</rng_seed>
        <r0>2.5</r0>
        <num_infected_seeds>540</num_infected_seeds>
        <seeding_age_min>1</seeding_age_min>
        <seeding_age_max>99</seeding_age_max>
        <start_date>2020-03-04</start_date>
        <stride_log_level>info</stride_log_level>
        <track_index_case>false</track_index_case>
        <use_install_dirs>true</use_install_dirs>
        <vaccine_profile>None</vaccine_profile>
        <num_daily_imported_cases>0</num_daily_imported_cases>
</run>
        )###";
}

ptree RunConfigManager::FromString(const string& s)
{
        ptree         pt;
        istringstream is(s);
        read_xml(is, pt, trim_whitespace);
        return pt;
}

string RunConfigManager::ToString(const ptree& pt)
{
        ostringstream ss;
        write_xml(ss, pt, xml_writer_make_settings<ptree::key_type>(' ', 8));
        return ss.str();
}

} // namespace util
} // namespace stride
