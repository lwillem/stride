
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
 *  Copyright 2020 Willem L, Kuylen E, Stijven S & Broeckhove J
 */

/**
 * @file
 * Implementation of scenario tests data.
 */

#include "ScenarioData.h"

#include "util/RunConfigManager.h"

#include <boost/property_tree/ptree.hpp>
#include <map>

using namespace std;
using namespace stride;
using namespace stride::util;
using boost::property_tree::ptree;

namespace Tests {

tuple<ptree, unsigned int, double> ScenarioData::Get(string tag)
{

	// Set model config (default = Measles)
	ptree pt = RunConfigManager::Create("TestsMeasles");


	// ... or set model config for influenza simulations
	if(tag.substr(0, 9) == "influenza"){
		pt = RunConfigManager::Create("TestsInfluenza");
	}

	// ... or set model config for covid19 simulations
	if(tag.substr(0, 7) == "covid19"){
		pt = RunConfigManager::Create("TestsCovid19");
	}

	// Set target values per scenario
	// note: last update on September 3rd, 2020, based on the mean of 20 realizations (MacOS, LW)
	const map<string, unsigned int> targets_default = {
		{"influenza_a", 500000U}, {"influenza_b", 0U}, {"influenza_c", 5U}, {"measles_16", 45000U},
		{"measles_26", 600000U},  {"r0_0", 1200U},     {"r0_4", 3400U},     {"r0_8", 9500U},
		{"r0_12", 23000U},        {"r0_16", 45000U},   {"covid19_base", 82500U}, {"covid19_all", 81000U},
		{"covid19_daily", 91000U},{"covid19_distancing", 19000U}, {"covid19_age_15min",90000U},
		{"covid19_householdclusters", 46000U}, {"covid19_tracing",41000U}, {"covid19_tracing_all",39000U},
		{"covid19_transm", 82500U},{"covid19_transm_gamma", 72300U},
		{"covid19_suscept", 82500U},{"covid19_suscept_age", 82500U},{"covid19_suscept_adapt", 56650U},
		{"covid19_fitting", 82500U},{"covid19_fitting_adapt", 41100U}};


	// Set margins per scenario
	const map<string, double> margins_default = {
		{"influenza_a", 1.0e-02}, {"influenza_b", 0.0}, {"influenza_c", 2.0e-02}, {"measles_16", 1.0e-01},
		{"measles_26", 5.0e-02},  {"r0_0", 5.0e-02},    {"r0_4", 1.0e-01},        {"r0_8", 1.0e-01},
		{"r0_12", 5.0e-02},       {"r0_16", 5.0e-02},   {"covid19_base", 1.0e-01},  {"covid19_all", 1.0e-01},
		{"covid19_daily", 1.0e-01},{"covid19_distancing", 1.0e-01},{"covid19_age_15min",1.0e-1},
		{"covid19_householdclusters", 1.5e-01}, // more stochastic effects observed
		{"covid19_tracing",1.0e-01}, {"covid19_tracing_all",1.0e-01},
		{"covid19_transm", 1.0e-01},{"covid19_transm_gamma", 1.0e-01},
		{"covid19_suscept", 1.0e-01},{"covid19_suscept_age", 1.0e-01},{"covid19_suscept_adapt", 1.0e-01},
		{"covid19_fitting", 1.0e-01},{"covid19_fitting_adapt", 1.0e-01}};

	unsigned int target;
	double       margin;
	target = targets_default.at(tag);
	margin = margins_default.at(tag);

	// Adjust some  parameters, according the scenario
	if (tag == "influenza_b") {
			pt.put("run.num_infected_seeds", 0);
	}
	if (tag == "influenza_c") {
			pt.put("run.num_infected_seeds", 5);
			pt.put("run.immunity_rate", 0.9991);
	}
	if (tag == "measles_16") {
			pt.put("run.r0", 16U);
	}
	if (tag == "measles_26") {
			pt.put("run.r0", 26U);
			pt.put("run.immunity_rate", 0U);
			pt.put("run.vaccine_rate", 0U);
			pt.put("run.num_days", 200U);
	}

	if (tag == "r0_0") {
			pt.put("run.r0", 0.0);
	}
	if (tag == "r0_4") {
			pt.put("run.r0", 4.0);
	}
	if (tag == "r0_8") {
			pt.put("run.r0", 8.0);
	}
	if (tag == "r0_12") {
			pt.put("run.r0", 12.0);
	}
	if (tag == "r0_16") {
			pt.put("run.r0", 16.0);
	}

	if (tag == "covid19_all") {
		pt.put("run.event_log_level", "All");
	}
	if (tag == "covid19_daily") {
			pt.put("run.num_daily_imported_cases", 10U);
	}
	if (tag == "covid19_distancing") {
			pt.put("run.holidays_file","calendar_gtester_covid19_distancing.csv");
			pt.put("run.cnt_reduction_workplace",1);
			pt.put("run.cnt_reduction_other",1);
//			pt.put("run.holidays_file","calendar_belgium_2020_covid19_exit_school_adjusted.csv");
//			pt.put("run.cnt_reduction_workplace",0.3);
//			pt.put("run.cnt_reduction_other",0.4);
//			pt.put("run.compliance_delay_workplace",2);
//			pt.put("run.compliance_delay_other",3);

	}
	if (tag == "covid19_age_15min") {
			pt.put("run.disease_config_file", "disease_covid19_age_15min.xml");
			pt.put("run.age_contact_matrix_file", "contact_matrix_flanders_conditional_teachers_15min.xml");
	}
	if (tag == "covid19_householdclusters") {
			pt.put("run.holidays_file", "calendar_belgium_2020_covid19_exit_schoolcategory_adjusted.csv");
			pt.put("run.start_date", "2020-06-01");
			pt.put("run.population_file", "pop_belgium600k_c500_teachers_censushh_extended3_size2.csv");
			pt.put("run.cnt_intensity_householdCluster", 4/7);
	}
	// set default tracing parameters
	if (tag == "covid19_tracing" || tag == "covid19_tracing_all") {
		    pt.put("run.event_log_level", "Transmissions");
			pt.put("run.holidays_file", "calendar_belgium_2020_covid19_exit_schoolcategory_adjusted.csv");
			pt.put("run.start_date", "2020-06-01");
			pt.put("run.detection_probability", 0.5);
			pt.put("run.tracing_efficiency_household", 1.0);
			pt.put("run.tracing_efficiency_other", 0.7);
			pt.put("run.test_false_negative", 0.1);
			pt.put("run.case_finding_capacity", 1000U);

	}
	// change log parameter to use the optimized version
	if (tag == "covid19_tracing_all") {
			pt.put("run.event_log_level", "ContactTracing");
	}

	if (tag == "covid19_transm") {
		pt.put("run.transmission_probability_distribution", "Constant");
		pt.put("run.transmission_probability_distribution_overdispersion", 0);
	}

	if (tag == "covid19_transm_gamma") {
		pt.put("run.transmission_probability_distribution", "Gamma");
		pt.put("run.transmission_probability_distribution_overdispersion", 0.8);
	}

	if (tag == "covid19_suscept") {
			pt.put("run.disease_susceptibility_age","1");
			pt.put("run.disease_susceptibility_agecat","0");
	}

	if (tag == "covid19_suscept_age") {
			pt.put("run.disease_susceptibility_age","1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,"
					"1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,"
					"1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1");
			pt.put("run.disease_susceptibility_agecat","0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,"
					"18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,"
					"44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,"
					"70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,"
					"96,97,98,99");
	}

	if (tag == "covid19_suscept_adapt") {
		pt.put("run.disease_susceptibility_age","1,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,1,0.9,0.9,"
				"0.9,0.9,0.9,0.9,0.9,0.9,1,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,1,0.9,0.9,0.9,0.9,"
				"0.9,0.9,0.9,0.9,1,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,1,0.9,0.9,0.9,0.9,0.9,0.9,"
				"0.9,0.9,1,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,1,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,"
				"1.0,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,1,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,1,0.9,"
				"0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9");
		pt.put("run.disease_susceptibility_agecat","0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,"
				"18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,"
				"44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,"
				"70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,"
				"96,97,98,99");
	}

	if (tag == "covid19_fitting") {
		pt.put("run.disease_susceptibility_age",0.05351622);
		pt.put("run.disease_susceptibility_agecat",0);
		pt.put("run.transmission_probability",1);
	}

	if (tag == "covid19_fitting_adapt") {
		pt.put("run.disease_susceptibility_age","0.02,0.053516219268665,0.08,0.053516219268665");
		pt.put("run.disease_susceptibility_agecat","0,18,59,70");
		pt.put("run.transmission_probability",1);
	}



	return make_tuple(pt, target, margin);
}

} // namespace Tests
