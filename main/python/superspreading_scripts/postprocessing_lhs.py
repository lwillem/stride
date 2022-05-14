############################################################################ #
#  This file is part of the Stride software.
#  It is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or any
#  later version.
#  The software is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  You should have received a copy of the GNU General Public License,
#  along with the software. If not, see <http://www.gnu.org/licenses/>.
#  see http://www.gnu.org/licenses/.
#
#
#  Copyright 2022, Kuylen E
############################################################################ #

import argparse
import matplotlib.pyplot as plt
import multiprocessing
import numpy as np
import os
import xml.etree.ElementTree as ET

from plots import plot_ar, plot_cumulative_cases_per_day, plot_day_of_last_infection, plot_extinction_probabilities, plot_final_size_frequencies, plot_new_cases_per_day, plot_p80s
from postprocessing_util import get_day_of_last_infection, get_experiment_ids, get_output_per_day, get_summary_output, get_total_cases

def main(output_dir, num_parallel_workers):
    scenario_names = ["baseline", "A", "B", "C", "D", "E", "F"]
    display_scenario_names = scenario_names

    extinction_threshold = 20

    num_days = np.nan
    population_size = np.nan

    all_days_last_infection = []
    all_days_last_infection_exclude_extinction = []

    all_final_sizes = []

    all_p80s = []
    all_p80s_exclude_extinction = []

    scenario_descriptions = {}

    for scenario_name in scenario_names:
        print(scenario_name)

        experiment_ids = get_experiment_ids(output_dir, scenario_name)
        with multiprocessing.Pool(processes=num_parallel_workers) as pool:
            config = ET.parse(os.path.join(output_dir, scenario_name, "exp{:04}".format(experiment_ids[0]), "config.xml")).getroot()
            alpha_i = float(config.find('transmission_probability_distribution_overdispersion').text)
            alpha_c = float(config.find('contact_distribution_overdispersion').text)

            scenario_descriptions[scenario_name] = {"alpha_i": alpha_i, "alpha_c": alpha_c}

            output_per_day = pool.starmap(get_output_per_day, [(output_dir, scenario_name, exp_id) for exp_id in experiment_ids])
            output_summary = pool.starmap(get_summary_output, [(output_dir, scenario_name, exp_id) for exp_id in experiment_ids])

            num_days = output_summary[0]["num_days"]
            population_size = output_summary[0]["population_size"]

            # Sort total cases from high to low & print
            # Used to determine extinction threshold
            total_cases = [get_total_cases(run["cases_per_day"], num_days) for run in output_per_day]
            total_cases.sort(reverse=True)
            print(total_cases)

            all_days_last_infection.append([get_day_of_last_infection(run["cases_per_day"], num_days) for run in output_per_day])
            all_days_last_infection_exclude_extinction.append([get_day_of_last_infection(run["cases_per_day"], num_days) for run in output_per_day if sum(run["cases_per_day"].values()) >= extinction_threshold])

            all_final_sizes.append(total_cases)

            # Calculate the proportion of infected individuals responsible for 80% of infections
            all_p80s.append([run["p80"] for run in output_summary])
            all_p80s_exclude_extinction.append([run["p80"] for run in output_summary if run["total_cases"] >= extinction_threshold])

            #plot_cumulative_cases_per_day(output_dir, "cumulative_cases_per_day", scenario_name, [run["cases_per_day"] for run in output_per_day], num_days, y_max=population_size)
            #plot_new_cases_per_day(output_dir, "new_cases_per_day", scenario_name, [run["cases_per_day"] for run in output_per_day], num_days, y_max=670000)

    display_scenario_names = ["baseline"] + ["ai = {:.2f},\n ac = {:.2f}".format(s["alpha_i"], s["alpha_c"]) for s in list(scenario_descriptions.values())[1:]]

    # Order by P80
    mean_p80s_exclude_extinction = [np.mean(p80s) for  p80s in all_p80s_exclude_extinction]
    indices_sorted_by_p80 = np.argsort(mean_p80s_exclude_extinction)[::-1]

    display_scenario_names_sorted = [display_scenario_names[i] for i in indices_sorted_by_p80]
    all_p80s_exclude_extinction_sorted = [all_p80s_exclude_extinction[i] for i in indices_sorted_by_p80]

    #plot_p80s(output_dir, "p80s", display_scenario_names_sorted, all_p80s, "blue")
    plot_p80s(output_dir, "p80s_exclude_extinction", display_scenario_names_sorted, all_p80s_exclude_extinction_sorted, "blue")

    #plot_ar(output_dir, "ar", display_scenario_names, all_final_sizes, num_days, population_size, violin_plot=True)
    #plot_ar(output_dir, "ar_exclude_extinction", display_scenario_names, all_final_sizes, num_days, population_size, extinction_threshold=extinction_threshold, y_min=0.8, y_max=1, violin_plot=True)

    #plot_day_of_last_infection(output_dir, "day_of_last_infection", display_scenario_names, all_days_last_infection, -1, num_days, violin_plot=True)
    #plot_day_of_last_infection(output_dir, "day_of_last_infection_exclude_extinction", display_scenario_names, all_days_last_infection_exclude_extinction, 100, num_days, violin_plot=True)

    #plot_final_size_frequencies(output_dir, "final_size_frequencies", display_scenario_names, all_final_sizes, "Outbreak size after {} days".format(num_days))
    #plot_extinction_probabilities(output_dir, "extinction_probabilities", display_scenario_names, all_final_sizes, extinction_threshold)


    # TODO peak size
    # TODO peak timing
    # TODO rt
    # TODO herd immunity threshold

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", type=str, help="Directory containing simulation results")
    parser.add_argument("--num_parallel_workers", type=int, default=4)

    args = parser.parse_args()
    main(args.output_dir, args.num_parallel_workers)


"""
from plots import plot_day_of_last_infection, plot_day_of_peak, plot_effective_r_by_day, plot_herd_immunity_threshold, plot_peak_sizes, plot_secondary_cases_distribution, plot_secondary_cases_per_index_case_histogram, plot_secondary_cases_per_index_case_means, plot_transmissions_by_location
from postprocessing_util import get_day_of_last_infection, get_herd_immunity_threshold, get_num_secondary_cases_frequencies, get_num_secondary_cases_per_index_case, get_transmissions_by_location

def main(output_dir, output_prefix, num_parallel_workers):

    for overdispersion_scenario in overdispersion_scenario_names:
        alpha = r"$\alpha_{i}$"
        if overdispersion_scenario == "contacts_overdispersion":
            alpha = r"$\alpha_{c}$"
        display_scenario_names = ["Baseline", alpha + " = 10", alpha + " = 1", alpha + " = 0.6", alpha + " = 0.4", alpha + " = 0.2"]

        all_cases_per_day = []

        all_hits = []
        all_hits_exclude_extinction = []
        all_hits_day = []
        all_hits_day_exclude_extinction = []

        all_secondary_cases = []
        all_secondary_cases_per_index_case = []

        for scenario_name in scenario_names:
            with multiprocessing.Pool(processes=num_parallel_workers) as pool:
                output_num_secondary_cases = pool.starmap(get_num_secondary_cases_frequencies, [(output_dir, scenario_name, exp_id) for exp_id in experiment_ids])
                secondary_cases_per_index_case = pool.starmap(get_num_secondary_cases_per_index_case, [(output_dir, scenario_name, exp_id) for exp_id in experiment_ids])
                transmissions_by_location = pool.starmap(get_transmissions_by_location, [(output_dir, scenario_name, exp_id) for exp_id in experiment_ids])


                all_cases_per_day.append([run["cases_per_day"] for run in output_per_day])

                all_secondary_cases_per_index_case.append(secondary_cases_per_index_case)


                all_hits.append([get_herd_immunity_threshold(run["rt_by_day"], run["cases_per_day"], num_days, population_size) for run in output_per_day])
                all_hits_exclude_extinction.append([get_herd_immunity_threshold(run["rt_by_day"], run["cases_per_day"], num_days, population_size) for run in output_per_day if sum(run["cases_per_day"].values()) >= extinction_threshold])

                all_hits_day.append([get_herd_immunity_threshold(run["rt_by_day"], run["cases_per_day"], num_days, population_size, get_day=True) for run in output_per_day])
                all_hits_day_exclude_extinction.append([get_herd_immunity_threshold(run["rt_by_day"], run["cases_per_day"], num_days, population_size, get_day=True) for run in output_per_day if sum(run["cases_per_day"].values()) >= extinction_threshold])

                all_secondary_cases.append(output_num_secondary_cases)

                plot_effective_r_by_day(output_dir, "rt", scenario_name, [run["rt_by_day"] for run in output_per_day], num_days, y_max=8, smoothed=True)


                plot_transmissions_by_location(output_dir, "transmissions_by_location", scenario_name, transmissions_by_location)

        plot_day_of_peak(output_dir, "day_of_peak_" + overdispersion_scenario, display_scenario_names, all_cases_per_day, 150, violin_plot=True)
        plot_day_of_peak(output_dir, "day_of_peak_exclude_extinction_" + overdispersion_scenario, display_scenario_names, all_cases_per_day, 150, extinction_threshold=extinction_threshold, violin_plot=True)

        plot_herd_immunity_threshold(output_dir, "hit_" + overdispersion_scenario, display_scenario_names, all_hits, show_day=False, violin_plot=True)
        plot_herd_immunity_threshold(output_dir, "hit_exclude_extinction_" + overdispersion_scenario, display_scenario_names, all_hits_exclude_extinction, show_day=False, y_min=0.4, violin_plot=True)

        plot_herd_immunity_threshold(output_dir, "hits_day_" + overdispersion_scenario, display_scenario_names, all_hits_day, show_day=True, violin_plot=True)
        plot_herd_immunity_threshold(output_dir, "hits_day_exclude_extinction_" + overdispersion_scenario, display_scenario_names, all_hits_day_exclude_extinction, show_day=True, violin_plot=True)

        plot_peak_sizes(output_dir, "peak_sizes_" + overdispersion_scenario, display_scenario_names, all_cases_per_day, 0, num_days, ymin=-10, ymax=670000, violin_plot=True)
        plot_peak_sizes(output_dir, "peak_sizes_exclude_extinction_" + overdispersion_scenario, display_scenario_names, all_cases_per_day, 0, num_days, extinction_threshold=extinction_threshold, ymin=400000, ymax=670000, violin_plot=True)

        plot_secondary_cases_distribution(output_dir, "secondary_cases_distribution_" + overdispersion_scenario, display_scenario_names, all_secondary_cases)

        plot_secondary_cases_per_index_case_histogram(output_dir, "secondary_cases_per_ic_" + overdispersion_scenario, all_secondary_cases_per_index_case, display_scenario_names)
        plot_secondary_cases_per_index_case_means(output_dir, "secondary_cases_per_ic_means_" + overdispersion_scenario, all_secondary_cases_per_index_case, display_scenario_names)

"""
