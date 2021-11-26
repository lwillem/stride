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
#  Copyright 2020, Willem L, Kuylen E & Broeckhove J
############################################################################ #

"""
"""

'''
import matplotlib.pyplot as plt
import numpy as np
import os

from util import get_secondary_cases_by_individual
from util import plot_ar, plot_day_last_infection, plot_final_size_frequencies

def plot_num_cases_over_period(output_dir, fig_name, all_num_cases, display_scenario_names, y_label):
    plt.violinplot(all_num_cases)
    plt.plot(range(1, len(all_num_cases) + 1), [np.mean(scen) for scen in all_num_cases], marker="o", linestyle="None")
    plt.xticks(range(1, len(all_num_cases) + 1), scenario_display_names, rotation=25)
    plt.ylabel(y_label)
    save_figure(output_dir, fig_name, extension="png")

def plot_resurgence_probabilities(output_dir, fig_name, all_num_case_after_release, display_scenario_names, extinction_threshold):
    resurgence_probs = []
    for scenario in all_num_case_after_release:
        resurgence_probs.append(len([run for run in scenario if run >= extinction_threshold]) / len(scenario))

    plt.bar(range(len(resurgence_probs)), resurgence_probs)
    plt.xticks(range(len(resurgence_probs)), display_scenario_names, rotation=25)
    plt.ylabel("Resurgence probability")
    plt.ylim((0, 1))
    save_figure(output_dir, fig_name)

def main(output_dir, scenario_names, display_scenario_names):
    all_total_cases = []

    all_last_days_with_infections = []

    for scenario in scenario_names:
        with multiprocessing.Pool(processes=4) as pool:

            # Calculate the final total number of cases for each simulation
            # Note: index cases are not counted.
            total_cases = [sum(x.values()) for x in secondary_cases_by_individual]
            all_total_cases.append(total_cases)


    plot_day_last_infection(output_dir, "day_of_last_infection", display_scenario_names, all_last_days_with_infections, violin_plot=True)
    plot_ar(output_dir, "ar", display_scenario_names, all_total_cases, num_days, pop_size, 0, violin_plot=True)
    plot_final_size_frequencies(output_dir, "final_size_frequencies", display_scenario_names, all_total_cases, num_days)

    plot_resurgence_probabilities(output_dir, "resurgence_probabilities", all_num_cases_release_exclude_extinct_before_d30, display_scenario_names, 500)

'''

"""import argparse
import multiprocessing

from plots import plot_cumulative_cases_per_day, plot_new_cases_per_day
from util import get_experiment_ids, get_cases_output, get_num_cases_over_period

def main(output_dir, scenario_names, display_scenario_names):
    num_days = 200

    if len(display_scenario_names) < len(scenario_names):
        display_scenario_names = scenario_names

    for s_i in range(len(scenario_names)):
        scenario_name = scenario_names[s_i]
        print(scenario_name)

        experiment_ids = get_experiment_ids(output_dir, scenario_name)
        with multiprocessing.Pool(processes=4) as pool:
            cases_output = pool.starmap(get_cases_output, [(output_dir, scenario_name, exp_id) for exp_id in experiment_ids])

            num_cases_before_lockdown = [get_num_cases_over_period(output["cases_by_day"], 0, 40) for output in cases_output]
            # TODO num cases during lockdown?
            # TODO num cases during release?

            plot_new_cases_per_day(output_dir, "sd_new_cases_per_day", scenario_name, [output["cases_by_day"] for output in cases_output], num_days)
            plot_cumulative_cases_per_day(output_dir, "sd_cumulative_cases_per_day", scenario_name, [output["cases_by_day"] for output in cases_output], num_days)

        # TODO plot num cases before / during / after lockdown?
'''
with multiprocessing.Pool(processes=4) as pool:
    # Calculate the proportion of infected individuals responsible for 80% of infections
    all_p80s.append([get_p80(output["secondary_cases_by_individual"], extinction_threshold) for output in cases_output]) # TODO exclude extinction?
    # Get final outbreak size
    total_cases = [get_total_cases(output["cases_by_day"], num_days) for output in cases_output]
    all_final_sizes.append(total_cases)
    # Day of last infection
    all_days_last_infection.append([get_day_of_last_infection(output["cases_by_day"], num_days) for output in cases_output])
    # TODO smoothed Rt by day?
    # Herd immunity threshold
    all_hits.append([get_herd_immunity_threshold(output["rt_by_day"], output["cases_by_day"], num_days, output["parameters"]["population_size"]) for output in cases_output])
    all_hits_day.append([get_herd_immunity_threshold(output["rt_by_day"], output["cases_by_day"], num_days, output["parameters"]["population_size"]) for output in cases_output], get_day=True)

    # Sort total cases from high to low & print
    # Used to determine extinction threshold
    total_cases.sort(reverse=True)
    print(total_cases)

    plot_effective_r_by_day(output_dir, "rt_by_day", scenario_name, [output["rt_by_day"] for output in cases_output], num_days)

plot_extinction_probabilities(output_dir, "extinction_probabilities", display_scenario_names, all_final_sizes, extinction_threshold)
plot_p80s(output_dir, "p80s", display_scenario_names, all_p80s)
plot_ar(output_dir, "ar", display_scenario_names, all_final_sizes, num_days, population_size)
plot_ar(output_dir, "ar_exclude_extinction", display_scenario_names, all_final_sizes, num_days, population_size, extinction_threshold=extinction_threshold)
plot_herd_immunity_threshold(output_dir, "hits", display_scenario_names, all_hits)
plot_herd_immunity_threshold(output_dir, "hits_day", display_scenario_names, all_hits_day, show_day=True)
plot_final_size_frequencies(output_dir, "final_size_frequencies", display_scenario_names, all_final_sizes, num_days)
plot_day_last_infection(output_dir, "last_day_with_infections", display_scenario_names, all_days_last_infection, num_days)
'''


import multiprocessing

from plots import plot_ar, plot_cumulative_cases_per_day, plot_day_of_last_infection, plot_day_of_peak, plot_extinction_probabilities, plot_effective_r_by_day, plot_final_size_frequencies, plot_herd_immunity_threshold, plot_new_cases_per_day, plot_p80s, plot_peak_sizes, plot_secondary_cases_distribution, plot_transmissions_by_location

from util import get_day_of_last_infection, get_experiment_ids, get_herd_immunity_threshold, get_output, get_p80, get_total_cases

def main(output_dir):
    for overdispersion_scenario in overdispersion_scenario_names:

        all_days_last_infection = []
        all_days_last_infection_exclude_extinction = []


        all_hits = []
        all_hits_exclude_extinction = []
        all_hits_day = []
        all_hits_day_exclude_extinction = []

        all_p80s = []
        all_p80s_exclude_extinction = []

        all_secondary_cases = []


        for scenario_name in scenario_names:

            with multiprocessing.Pool(processes=4) as pool:


                all_days_last_infection.append([get_day_of_last_infection(run["cases_per_day"], num_days) for run in output])
                all_days_last_infection_exclude_extinction.append([get_day_of_last_infection(run["cases_per_day"], num_days) for run in output if sum(run["cases_per_day"].values()) >= extinction_threshold])

                all_secondary_cases.append([run["secondary_cases_by_individual"] for run in output])

                # Calculate the proportion of infected individuals responsible for 80% of infections
                all_p80s.append([get_p80(run["secondary_cases_by_individual"]) for run in output])
                all_p80s_exclude_extinction.append([get_p80(run["secondary_cases_by_individual"], extinction_threshold = extinction_threshold) for run in output])

                all_hits.append([get_herd_immunity_threshold(run["rt_by_day"], run["cases_per_day"], num_days, run["parameters"]["population_size"]) for run in output])
                all_hits_exclude_extinction.append([get_herd_immunity_threshold(run["rt_by_day"], run["cases_per_day"], num_days, run["parameters"]["population_size"]) for run in output if sum(run["cases_per_day"].values()) >= extinction_threshold])

                all_hits_day.append([get_herd_immunity_threshold(run["rt_by_day"], run["cases_per_day"], num_days, run["parameters"]["population_size"], get_day=True) for run in output])
                all_hits_day_exclude_extinction.append([get_herd_immunity_threshold(run["rt_by_day"], run["cases_per_day"], num_days, run["parameters"]["population_size"], get_day=True) for run in output if sum(run["cases_per_day"].values()) >= extinction_threshold])

                plot_effective_r_by_day(output_dir, "rt", scenario_name, [run["rt_by_day"] for run in output], num_days, y_max=30)

                plot_transmissions_by_location(output_dir, "transmissions_by_location", scenario_name, [run["transmissions_by_location"] for run in output])




        plot_day_of_peak(output_dir, "day_of_peak_" + overdispersion_scenario, display_scenario_names, all_cases_per_day, num_days)
        plot_day_of_peak(output_dir, "day_of_peak_exclude_extinction_" + overdispersion_scenario, display_scenario_names, all_cases_per_day, num_days, extinction_threshold=extinction_threshold)

        plot_day_of_last_infection(output_dir, "day_of_last_infection_" + overdispersion_scenario, display_scenario_names, all_days_last_infection, num_days)
        plot_day_of_last_infection(output_dir, "day_of_last_infection_exclude_extinction_" + overdispersion_scenario, display_scenario_names, all_days_last_infection_exclude_extinction, num_days)

        plot_herd_immunity_threshold(output_dir, "hit_" + overdispersion_scenario, display_scenario_names, all_hits, show_day=False)
        plot_herd_immunity_threshold(output_dir, "hits_day_" + overdispersion_scenario, display_scenario_names, all_hits_day, show_day=True)

        plot_herd_immunity_threshold(output_dir, "hit_exclude_extinction_" + overdispersion_scenario, display_scenario_names, all_hits_exclude_extinction, show_day=False)
        plot_herd_immunity_threshold(output_dir, "hits_day_exclude_extinction_" + overdispersion_scenario, display_scenario_names, all_hits_day_exclude_extinction, show_day=True)

        plot_final_size_frequencies(output_dir, "final_size_frequencies_" + overdispersion_scenario, display_scenario_names, all_final_sizes, num_days)
        plot_extinction_probabilities(output_dir, "extinction_probabilities_" + overdispersion_scenario, display_scenario_names, all_final_sizes, extinction_threshold)

        plot_p80s(output_dir, "p80s_" + overdispersion_scenario, display_scenario_names, all_p80s)
        plot_p80s(output_dir, "p80s_exclude_extinction_" + overdispersion_scenario, display_scenario_names, all_p80s_exclude_extinction)

        plot_peak_sizes(output_dir, "peak_sizes_" + overdispersion_scenario, display_scenario_names, all_cases_per_day, ymin=0, ymax=175000)
        plot_peak_sizes(output_dir, "peak_sizes_exclude_extinction_" + overdispersion_scenario, display_scenario_names, all_cases_per_day, extinction_threshold=extinction_threshold, ymin=100000, ymax=175000)

        plot_secondary_cases_distribution(output_dir, "secondary_cases_distribution_" + overdispersion_scenario, display_scenario_names, all_secondary_cases)

"""

import argparse
import multiprocessing
import numpy as np

from plots import plot_ar, plot_cumulative_cases_per_day, plot_new_cases_per_day, plot_num_cases_over_period
from util import get_experiment_ids, get_output, get_total_cases

def main(output_dir):
    baseline_scenario_name = "sd_baseline"

    overdispersion_scenario_names = ["sd_infectiousness_overdispersion", "sd_contacts_overdispersion"]
    overdispersion_parameters = ["1000", "100", "60", "40", "20"]

    num_days = 200
    population_size = np.nan

    for overdispersion_scenario in overdispersion_scenario_names:
        scenario_names = [baseline_scenario_name] + [overdispersion_scenario + "_" + overdispersion for overdispersion in overdispersion_parameters]

        alpha = r"$\alpha$"
        display_scenario_names = ["Baseline", alpha + " = 10", alpha + " = 1", alpha + " = 0.6", alpha + " = 0.4", alpha + " = 0.2"]

        all_cases_per_day = []
        all_final_sizes = []

        for scenario_name in scenario_names:
            print(scenario_name)

            experiment_ids = get_experiment_ids(output_dir, scenario_name)
            output = {}
            with multiprocessing.Pool(processes=4) as pool:
                output = pool.starmap(get_output, [(output_dir, scenario_name, exp_id) for exp_id in experiment_ids])

                all_cases_per_day.append([run["cases_per_day"] for run in output])

                # Sort total cases from high to low & print
                # Used to determine extinction threshold
                total_cases = [get_total_cases(run["cases_per_day"], num_days) for run in output]
                total_cases.sort(reverse=True)
                print(total_cases)

                all_final_sizes.append(total_cases)

                plot_new_cases_per_day(output_dir, "new_cases_per_day", scenario_name, [run["cases_per_day"] for run in output], num_days, y_max=170000)
                plot_cumulative_cases_per_day(output_dir, "cumulative_cases_per_day", scenario_name, [run["cases_per_day"] for run in output], num_days, y_max=3000000)

                population_size = output[0]["parameters"]["population_size"]


        #plot_ar(output_dir, "ar_" + overdispersion_scenario, display_scenario_names, all_final_sizes, num_days, population_size)
        plot_num_cases_over_period(output_dir, "cases_before_lockdown_" + overdispersion_scenario, display_scenario_names, 0, 30, all_cases_per_day)
        plot_num_cases_over_period(output_dir, "cases_during_lockdown_" + overdispersion_scenario, display_scenario_names, 30, 90, all_cases_per_day)
        plot_num_cases_over_period(output_dir, "cases_after_lockdown_" + overdispersion_scenario, display_scenario_names, 90, 120, all_cases_per_day)

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", type=str)

    args = parser.parse_args()
    main(args.output_dir)
