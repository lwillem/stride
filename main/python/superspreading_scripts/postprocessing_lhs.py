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

from postprocessing_util import get_day_of_last_infection, get_experiment_ids, get_herd_immunity_threshold, get_output_per_day, get_num_secondary_cases_per_index_case, get_summary_output, get_total_cases
from plots import plot_secondary_cases_per_index_case_means, save_figure

def get_peak_day(cases_per_day):
    peak_day = -1
    peak_cases = -1
    for day, num_cases in cases_per_day.items():
        if num_cases > peak_cases:
            peak_cases = num_cases
            peak_day = day
    return peak_day

def plot_contours(output_dir, figure_name, alpha_is, alpha_cs, means, min, max):
    alpha_i_str = r"$\alpha_{i}$"
    alpha_c_str = r"$\alpha_{c}$"

    levels = np.linspace(min, max, num=9)

    contourplot = plt.tricontourf(alpha_is, alpha_cs, means, cmap=cm.coolwarm, levels=levels)
    plt.plot(alpha_is, alpha_cs, 'o', markersize=2, color="grey")

    plt.xlabel(alpha_i_str)
    plt.ylabel(alpha_c_str)
    plt.xlim(0.2,0.6)
    plt.ylim(0.2,0.6)
    plt.colorbar(contourplot)

    save_figure(output_dir, figure_name)

def main(output_dir, num_scenarios, num_parallel_workers):
    num_days = np.nan
    population_size = np.nan

    extinction_threshold = 50

    alpha_is = []
    alpha_cs = []

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

    for s_i in range(num_scenarios):
        scenario_name = "scenario_" + str(s_i)
        experiment_ids = get_experiment_ids(output_dir, scenario_name)
        with multiprocessing.Pool(processes=num_parallel_workers) as pool:
            config = ET.parse(os.path.join(output_dir, scenario_name, "exp{:04}".format(experiment_ids[0]), "config.xml")).getroot()
            alpha_is.append(float(config.find('transmission_probability_distribution_overdispersion').text))
            alpha_cs.append(float(config.find('contact_distribution_overdispersion').text))

            output_per_day = pool.starmap(get_output_per_day, [(output_dir, scenario_name, exp_id) for exp_id in experiment_ids])
            secondary_cases_per_index_case = pool.starmap(get_num_secondary_cases_per_index_case, [(output_dir, scenario_name, exp_id) for exp_id in experiment_ids])
            output_summary = pool.starmap(get_summary_output, [(output_dir, scenario_name, exp_id) for exp_id in experiment_ids])

            num_days = output_summary[0]["num_days"]
            population_size = output_summary[0]["population_size"]

            all_cases_per_day.append([run["cases_per_day"] for run in output_per_day])
            all_days_last_infection.append([get_day_of_last_infection(run["cases_per_day"], num_days) for run in output_per_day])
            all_days_last_infection_exclude_extinction.append([get_day_of_last_infection(run["cases_per_day"], num_days) for run in output_per_day if sum(run["cases_per_day"].values()) >= extinction_threshold])

            # Sort total cases from high to low & print
            # Used to determine extinction threshold
            total_cases = [get_total_cases(run["cases_per_day"], num_days) for run in output_per_day]
            total_cases.sort(reverse=True)
            print(total_cases)

            # Calculate the proportion of infected individuals responsible for 80% of infections
            all_p80s.append([run["p80"] for run in output_summary])
            all_p80s_exclude_extinction.append([run["p80"] for run in output_summary if run["total_cases"] >= extinction_threshold])

            all_hits.append([get_herd_immunity_threshold(run["rt_by_day"], run["cases_per_day"], num_days, population_size) for run in output_per_day])
            all_hits_exclude_extinction.append([get_herd_immunity_threshold(run["rt_by_day"], run["cases_per_day"], num_days, population_size) for run in output_per_day if sum(run["cases_per_day"].values()) >= extinction_threshold])
            all_hits_day.append([get_herd_immunity_threshold(run["rt_by_day"], run["cases_per_day"], num_days, population_size, get_day=True) for run in output_per_day])
            all_hits_day_exclude_extinction.append([get_herd_immunity_threshold(run["rt_by_day"], run["cases_per_day"], num_days, population_size, get_day=True) for run in output_per_day if sum(run["cases_per_day"].values()) >= extinction_threshold])

            all_final_sizes.append(total_cases)

            all_secondary_cases_per_index_case.append(secondary_cases_per_index_case)

    # Plot P80s
    #mean_p80s = [np.nanmean(scenario) for scenario in all_p80s]
    mean_p80s_exclude_extinction = [np.nanmean(scenario) for scenario in all_p80s_exclude_extinction]
    #plot_contours(output_dir, "lhs_p80s", alpha_is, alpha_cs, mean_p80s, 0, 1)
    plot_contours(output_dir, "lhs_p80s_exclude_extinction", alpha_is, alpha_cs, mean_p80s_exclude_extinction, 0.05, 0.2)

    # Plot extinction probabilities
    extinction_probabilities = [sum(1 for c in scenario if c < extinction_threshold) / len(scenario) for scenario in all_final_sizes]
    plot_contours(output_dir, "lhs_extinction_probabilities", alpha_is, alpha_cs, extinction_probabilities, 0.4, 0.8)

    # Plot attack rate
    #mean_ars = [np.nanmean([c / population_size for c in scenario]) for scenario in all_final_sizes]
    mean_ars_exclude_extinction = [np.nanmean([c / population_size for c in scenario if c >= extinction_threshold]) for scenario in all_final_sizes]
    #plot_contours(output_dir, "lhs_ar", alpha_is, alpha_cs, mean_ars, 0, 1)
    plot_contours(output_dir, "lhs_ar_exclude_extinction", alpha_is, alpha_cs, mean_ars_exclude_extinction, 0.8, 0.9)

    # Plot herd immunity threshold
    #mean_hits = [np.nanmean(scenario) for scenario in all_hits]
    mean_hits_exclude_extinction = [np.nanmean(scenario) for scenario in all_hits_exclude_extinction]
    #plot_contours(output_dir, "lhs_hit", alpha_is, alpha_cs, mean_hits, 0, 1)
    plot_contours(output_dir, "lhs_hit_exclude_extinction", alpha_is, alpha_cs, mean_hits_exclude_extinction, 0.55, 0.65)

    # Plot peak sizes
    #mean_peak_sizes = [np.nanmean([max(c.values()) for c in scenario]) for scenario in all_cases_per_day]
    mean_peak_sizes_exclude_extinction = [np.nanmean([max(c.values()) for c in scenario if sum(c.values()) >= extinction_threshold]) for scenario in all_cases_per_day]
    #plot_contours(output_dir, "lhs_peak_sizes", alpha_is, alpha_cs, mean_peak_sizes, 0, 600000)
    plot_contours(output_dir, "lhs_peak_sizes_exclude_extinction", alpha_is, alpha_cs, mean_peak_sizes_exclude_extinction, 460000, 540000)

    # Plot peak timing
    #mean_day_of_peak = [np.nanmean([get_peak_day(run) for run in scenario]) for scenario in all_cases_per_day]
    mean_day_of_peak_exclude_extinction = [np.nanmean([get_peak_day(run) for run in scenario if sum(run.values()) >= extinction_threshold]) for scenario in all_cases_per_day]
    #plot_contours(output_dir, "lhs_day_of_peak", alpha_is, alpha_cs, mean_day_of_peak, 0, 200)
    plot_contours(output_dir, "lhs_day_of_peak_exclude_extinction", alpha_is, alpha_cs, mean_day_of_peak_exclude_extinction, 60, 80)

    # Plot day on which herd immunity threshold is reached
    #mean_hits_day = [np.nanmean(scenario) for scenario in all_hits_day]
    mean_hits_day_exclude_extinction = [np.nanmean(scenario) for scenario in all_hits_day_exclude_extinction]
    #plot_contours(output_dir, "lhs_hit_day", alpha_is, alpha_cs, mean_hits_day, 0, 200)
    plot_contours(output_dir, "lhs_hit_day_exclude_extinction", alpha_is, alpha_cs, mean_hits_day_exclude_extinction, 60, 90)

    # Plot day of last infection
    #mean_day_last_infection = [np.nanmean(scenario) for scenario in all_days_last_infection]
    mean_day_last_infection_exclude_extinction = [np.nanmean(scenario) for scenario in all_days_last_infection_exclude_extinction]
    #plot_contours(output_dir, "lhs_day_of_last_infection", alpha_is, alpha_cs, mean_day_last_infection, 0, 200)
    plot_contours(output_dir, "lhs_day_of_last_infection_exclude_extinction", alpha_is, alpha_cs, mean_day_last_infection_exclude_extinction, 140, 160)

    plot_secondary_cases_per_index_case_means(output_dir, "lhs_secondary_cases_per_ic", all_secondary_cases_per_index_case, list(range(num_scenarios)), "blue", x_label="Scenario number", label_rotation=90, xtick_step=3)

if __name__=="__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("output_dir", type=str, help="Directory containing simulation results.")
    parser.add_argument("num_scenarios", type=int, help="Number of combinations in LHS grid.")
    parser.add_argument("--num_parallel_workers", type=int, default=4)

    args = parser.parse_args()
    main(args.output_dir, args.num_scenarios, args.num_parallel_workers)
