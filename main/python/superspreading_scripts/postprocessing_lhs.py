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

from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

from plots import plot_ar, plot_day_of_last_infection, plot_day_of_peak, plot_effective_r_by_day, plot_extinction_probabilities, plot_final_size_frequencies, plot_herd_immunity_threshold, plot_new_cases_per_day, plot_p80s, plot_peak_sizes, plot_secondary_cases_per_index_case_means, save_figure
from postprocessing_util import get_day_of_last_infection, get_experiment_ids, get_herd_immunity_threshold, get_num_secondary_cases_per_index_case, get_output_per_day, get_summary_output, get_total_cases

def plot_contours(alpha_is, alpha_cs, means):
    alpha_i_str = r"$\alpha_{i}$"
    alpha_c_str = r"$\alpha_{c}$"

    contourplot = plt.tricontourf(alpha_is, alpha_cs, means, cmap=cm.coolwarm)
    plt.plot(alpha_is, alpha_cs, 'o', markersize=2, color="grey")
    plt.colorbar(contourplot)
    plt.xlabel(alpha_i_str)
    plt.ylabel(alpha_c_str)

    plt.show()

def main(output_dir, num_parallel_workers):
    display_scenario_names = scenario_names

    extinction_threshold = 20

    num_days = np.nan
    population_size = np.nan

    all_cases_per_day = []

    all_days_last_infection = []
    all_days_last_infection_exclude_extinction = []

    all_final_sizes = []

    all_hits = []
    all_hits_exclude_extinction = []
    all_hits_day = []
    all_hits_day_exclude_extinction = []

    all_p80s = []
    all_p80s_exclude_extinction = []

    all_secondary_cases_per_index_case = []

    #scenario_descriptions = {}

    xs = []
    ys = []

    for scenario_name in scenario_names:
        print(scenario_name)

        experiment_ids = get_experiment_ids(output_dir, scenario_name)
        with multiprocessing.Pool(processes=num_parallel_workers) as pool:
            config = ET.parse(os.path.join(output_dir, scenario_name, "exp{:04}".format(experiment_ids[0]), "config.xml")).getroot()
            alpha_i = float(config.find('transmission_probability_distribution_overdispersion').text)
            alpha_c = float(config.find('contact_distribution_overdispersion').text)
            xs.append(alpha_i)
            ys.append(alpha_c)

            #scenario_descriptions[scenario_name] = {"alpha_i": alpha_i, "alpha_c": alpha_c}

            output_per_day = pool.starmap(get_output_per_day, [(output_dir, scenario_name, exp_id) for exp_id in experiment_ids])
            output_summary = pool.starmap(get_summary_output, [(output_dir, scenario_name, exp_id) for exp_id in experiment_ids])

            secondary_cases_per_index_case = pool.starmap(get_num_secondary_cases_per_index_case, [(output_dir, scenario_name, exp_id) for exp_id in experiment_ids])

            num_days = output_summary[0]["num_days"]
            population_size = output_summary[0]["population_size"]

            all_cases_per_day.append([run["cases_per_day"] for run in output_per_day])
            all_secondary_cases_per_index_case.append(secondary_cases_per_index_case)

            # Sort total cases from high to low & print
            # Used to determine extinction threshold
            total_cases = [get_total_cases(run["cases_per_day"], num_days) for run in output_per_day]
            total_cases.sort(reverse=True)
            print(total_cases)

            all_final_sizes.append(total_cases)

            #all_days_last_infection.append([get_day_of_last_infection(run["cases_per_day"], num_days) for run in output_per_day])
            #all_days_last_infection_exclude_extinction.append([get_day_of_last_infection(run["cases_per_day"], num_days) for run in output_per_day if sum(run["cases_per_day"].values()) >= extinction_threshold])

            #all_hits.append([get_herd_immunity_threshold(run["rt_by_day"], run["cases_per_day"], num_days, population_size) for run in output_per_day])
            #all_hits_exclude_extinction.append([get_herd_immunity_threshold(run["rt_by_day"], run["cases_per_day"], num_days, population_size) for run in output_per_day if sum(run["cases_per_day"].values()) >= extinction_threshold])

            #all_hits_day.append([get_herd_immunity_threshold(run["rt_by_day"], run["cases_per_day"], num_days, population_size, get_day=True) for run in output_per_day])
            #all_hits_day_exclude_extinction.append([get_herd_immunity_threshold(run["rt_by_day"], run["cases_per_day"], num_days, population_size, get_day=True) for run in output_per_day if sum(run["cases_per_day"].values()) >= extinction_threshold])

            # Calculate the proportion of infected individuals responsible for 80% of infections
            all_p80s.append([run["p80"] for run in output_summary])
            #all_p80s_exclude_extinction.append([run["p80"] for run in output_summary if run["total_cases"] >= extinction_threshold])



    # Order by P80
    #mean_p80s_exclude_extinction = [np.mean(p80s) for  p80s in all_p80s_exclude_extinction]
    #indices_sorted_by_p80 = np.argsort(mean_p80s_exclude_extinction)[::-1]

    #display_scenario_names_sorted = [display_scenario_names[i] for i in indices_sorted_by_p80]

    #plot_p80s(output_dir, "p80s", display_scenario_names_sorted, all_p80s, "blue")
    #plot_p80s(output_dir, "p80s_exclude_extinction", display_scenario_names_sorted, [all_p80s_exclude_extinction[i] for i in indices_sorted_by_p80], "blue")

    #plot_ar(output_dir, "ar", display_scenario_names_sorted, [all_final_sizes[i] for i in indices_sorted_by_p80], num_days, population_size, "blue")
    #plot_ar(output_dir, "ar_exclude_extinction", display_scenario_names_sorted, [all_final_sizes[i] for i in indices_sorted_by_p80], num_days, population_size, "blue", extinction_threshold=extinction_threshold, y_min=0.8, y_max=1)

    mean_p80s = [np.nanmean(p80s) for p80s in all_p80s]
    #plot_contours(xs, ys, mean_p80s)

    mean_ars = [np.nanmean([c / population_size for c in final_sizes]) for final_sizes in all_final_sizes]
    print(mean_ars)

    '''mean_final_sizes = [np.mean(result) for result in all_final_sizes][1:]

    contourplot = plt.tricontourf(xs[1:], ys[1:], [np.mean(result) for result in all_final_sizes][1:], cmap=cm.coolwarm)
    plt.plot(xs[1:], ys[1:], 'o', markersize=2, color='grey')

    cbar = plt.colorbar(contourplot)
    plt.xlabel(alpha_i_str)
    plt.ylabel(alpha_c_str)
    plt.show()'''


    #plot_day_of_last_infection(output_dir, "day_of_last_infection", display_scenario_names_sorted, [all_days_last_infection[i] for i in indices_sorted_by_p80], -1, num_days, "blue")
    #plot_day_of_last_infection(output_dir, "day_of_last_infection_exclude_extinction", display_scenario_names_sorted, [all_days_last_infection_exclude_extinction[i] for i in indices_sorted_by_p80], 100, num_days, "blue")

    #plot_final_size_frequencies(output_dir, "final_size_frequencies", display_scenario_names_sorted, [all_final_sizes[i] for i in indices_sorted_by_p80], "Outbreak size after {} days".format(num_days))

    #plot_day_of_peak(output_dir, "day_of_peak", display_scenario_names_sorted, [all_cases_per_day[i] for i in indices_sorted_by_p80], 150, "blue")
    #plot_day_of_peak(output_dir, "day_of_peak_exclude_extinction", display_scenario_names_sorted, [all_cases_per_day[i] for i in indices_sorted_by_p80], 150, "blue", extinction_threshold=extinction_threshold)

    #plot_extinction_probabilities(output_dir, "extinction_probabilities", display_scenario_names_sorted, [all_final_sizes[i] for i in indices_sorted_by_p80], extinction_threshold, "blue")

    #plot_peak_sizes(output_dir, "peak_sizes", display_scenario_names_sorted, [all_cases_per_day[i] for i in indices_sorted_by_p80], 0, num_days, "blue", ymin=-10, ymax=670000)
    #plot_peak_sizes(output_dir, "peak_sizes_exclude_extinction", display_scenario_names_sorted, [all_cases_per_day[i] for i in indices_sorted_by_p80], 0, num_days, "blue", extinction_threshold=extinction_threshold, ymin=400000, ymax=670000)

    #plot_herd_immunity_threshold(output_dir, "hit", display_scenario_names_sorted, [all_hits[i] for i in indices_sorted_by_p80], "blue", show_day=False)
    #plot_herd_immunity_threshold(output_dir, "hit_exclude_extinction", display_scenario_names_sorted, [all_hits_exclude_extinction[i] for i in indices_sorted_by_p80], "blue", show_day=False, y_min=0.4)

    #plot_herd_immunity_threshold(output_dir, "hits_day", display_scenario_names_sorted, [all_hits_day[i] for i in indices_sorted_by_p80], "blue", show_day=True)
    #plot_herd_immunity_threshold(output_dir, "hits_day_exclude_extinction", display_scenario_names_sorted, [all_hits_day_exclude_extinction[i] for i in indices_sorted_by_p80], "blue", show_day=True)

    #plot_secondary_cases_per_index_case_means(output_dir, "secondary_cases_per_ic_means", [all_secondary_cases_per_index_case[i] for i in indices_sorted_by_p80], display_scenario_names_sorted, "blue")


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", type=str, help="Directory containing simulation results.")
    parser.add_argument("num_scenarios", type=int, help="Number of combinations in LHS grid.")
    parser.add_argument("--num_parallel_workers", type=int, default=4)

    args = parser.parse_args()
    main(args.output_dir, args.num_scenarios, args.num_parallel_workers)
