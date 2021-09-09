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
import multiprocessing
import numpy as np
import os

from statsmodels.nonparametric.smoothers_lowess import lowess

from util import get_cumulative_cases, get_new_cases_per_day, get_secondary_cases_by_individual, save_figure
from util import plot_ar, plot_day_last_infection, plot_evolution, plot_final_size_frequencies

def get_p80(secondary_cases_by_individual, total_cases, extinction_threshold=0):
    """
        Get proportion of infected individuals responsible for 80% of infections.
    """
    p80s = []
    for run_i in range(len(secondary_cases_by_individual)):
        total_cases_run = total_cases[run_i]
        if total_cases_run == 0: # If no transmissions occur, P80 cannot be calculated
            pass
        if total_cases_run >= extinction_threshold:
            # Order number of secondary cases by individual in decreasing order
            secondary_cases_sorted = list(secondary_cases_by_individual[run_i].values())
            secondary_cases_sorted.sort(reverse=True)

            num_cases_responsible = 0
            num_cases_caused = 0
            for s in secondary_cases_sorted:
                num_cases_caused += s
                num_cases_responsible += 1
                if num_cases_caused >= (total_cases_run * 0.80):
                    break
            p80s.append(num_cases_responsible / total_cases_run)
    return p80s

def plot_p80s(output_dir, fig_name, display_scenario_names, all_p80s):
    plt.boxplot(all_p80s, labels=display_scenario_names)
    plt.xticks(rotation=25)
    plt.ylabel("P80")

    save_figure(output_dir, fig_name)

def plot_extinction_probabilities(output_dir, fig_name, display_scenario_names, all_total_cases, extinction_threshold=0):
    extinction_probabilities = []
    for scenario_total_cases in all_total_cases:
        extinction_probabilities.append(len([x for x in scenario_total_cases if x < extinction_threshold]) / len(scenario_total_cases))

    plt.bar(range(len(all_total_cases)), extinction_probabilities)
    plt.xticks(range(len(display_scenario_names)), display_scenario_names, rotation=25)
    plt.ylabel("Extinction probability (threshold = {} cases)".format(extinction_threshold))
    plt.ylim(0, 1.1)

    save_figure(output_dir, fig_name)

def get_effective_r_by_day(output_dir, scenario_name, experiment_id, num_days, extinction_threshold=0, fraction=False):
    infected_by_day = {} # Keep track of IDs of persons infected per day
    potential_infectors = {} # Keep track of number of infections caused by infected indviduals

    for day in range(num_days):
        infected_by_day[day] = []

    log_file = os.path.join(output_dir, scenario_name, "exp" + "{:04}".format(experiment_id), "event_log.txt")
    with open(log_file) as f:
        for line in f:
            line = line.split(" ")
            tag = line[0]
            if tag == "[PRIM]": # Index case
                infected_id = int(float(line[1]))

                sim_day = int(line[6])
                infected_by_day[sim_day].append(infected_id)

                # Add infected to potential infectors
                if infected_id not in potential_infectors:
                    potential_infectors[infected_id] = 0

            elif tag == "[TRAN]": # Transmission
                infected_id = int(float(line[1]))
                sim_day = int(line[6])
                infected_by_day[sim_day].append(infected_id)

                # Add infected to potential infectors
                if infected_id not in potential_infectors:
                    potential_infectors[infected_id] = 0
                # Add infector to potential infectors
                # and add this transmission to total infections count
                infector_id = int(float(line[2]))
                if infector_id not in potential_infectors:
                    potential_infectors[infector_id] = 1
                else:
                    potential_infectors[infector_id] += 1

    # Only return result if not extinction case
    if sum(potential_infectors.values()) < extinction_threshold:
        return None

    rt_by_day = {}
    for day in infected_by_day:
        if len(infected_by_day[day]) == 0: # No individuals infected on this day -> no Rt
            rt_by_day[day] = np.nan
        else:
            # Get number of secondary cases caused by every individual infected on this day
            secondary_cases_on_day = []
            for infector_id in infected_by_day[day]:
                secondary_cases_on_day.append(potential_infectors[infector_id])

            if fraction: # Get fraction of Rts >= 1 on day
                rt_by_day[day] = len([x for x in secondary_cases_on_day if x >= 1]) / len(secondary_cases_on_day)
            else: # Get mean Rt on day
                rt_by_day[day] = sum(secondary_cases_on_day) / len(secondary_cases_on_day)

    return rt_by_day

def get_herd_immunity_threshold(output_dir, scenario_name, experiment_id, num_days, pop_size, extinction_threshold=0, get_day=False):
    rt_by_day = get_effective_r_by_day(output_dir, scenario_name, experiment_id, num_days, extinction_threshold)

    if rt_by_day is None: # Do not calculate number for extinction cases
        return np.nan

    # Smooth (using LOWESS function)
    rt_by_day_smoothed = lowess([rt_by_day[day] for day in range(num_days)], range(num_days), is_sorted=True, return_sorted=False)

    # Look for last day where Rt >= 1
    last_day_rt_greq_1 = np.nan

    for day in range(num_days - 1, -1, -1): # Iterate over simulation days in reverse order
        #rt = rt_by_day[day]
        rt = rt_by_day_smoothed[day]
        if rt >= 1:
            last_day_rt_greq_1 = day
            break

    if get_day: # Return day on which herd immunity threshold is reached
        return last_day_rt_greq_1
    else:
        # Get proportion of population no longer susceptible on this day
        # = cumulative_cases (including index cases) / population size
        cumulative_cases = get_cumulative_cases(output_dir, scenario_name, experiment_id, last_day_rt_greq_1 + 1, include_index_cases=True)
        herd_immunity_threshold = cumulative_cases / pop_size

        return herd_immunity_threshold

def plot_effective_r_by_day(output_dir, fig_name, scenario_name, rt_by_day, num_days, fraction=False):
    mean = []
    lower = []
    upper = []

    for day in range(num_days):
        er_on_day = [run[day] for run in rt_by_day]
        er_on_day = [result for result in er_on_day if not np.isnan(result)]

        mean.append(np.mean(er_on_day) if len(er_on_day) > 0 else np.nan)
        lower.append(np.percentile(er_on_day, 2.5) if len(er_on_day) > 0 else np.nan)
        upper.append(np.percentile(er_on_day, 97.5) if len(er_on_day) > 0 else np.nan)

    plt.plot(range(num_days), mean)
    plt.fill_between(range(num_days), lower, upper, color="lightgrey")

    if not fraction: # If looking at mean Rt per day, add reference line at Rt = 1
        plt.plot(range(num_days), [1] * num_days)

    plt.xlabel("Simulation day")
    plt.xlim(0, num_days)

    if fraction:
        plt.ylabel("Fraction of Rt >= 1")
        plt.ylim(-0.1, 1.1)
    else:
        plt.ylabel("Mean Rt")
        plt.ylim(-0.5, 22)

    save_figure(output_dir, fig_name + "_" + scenario_name)

def plot_herd_immunity_threshold(output_dir, fig_name, display_scenario_names, all_herd_immunity_thresholds, show_day=False):
    means = [np.mean(scenario_results) for scenario_results in all_herd_immunity_thresholds]

    plt.violinplot(all_herd_immunity_thresholds)
    plt.plot(range(1, len(means) + 1), means, linestyle="None", marker="o")
    plt.xticks(range(1, len(scenario_display_names) + 1), display_scenario_names, rotation=25)

    if show_day:
        plt.ylabel("Day on which Rt >= 1 for the last time")
    else:
        plt.ylabel("Herd immunity threshold")

    save_figure(output_dir, fig_name, extension="png") # Use png extension to accurately save opaqueness in figure

def main(output_dir, scenario_names, display_scenario_names):
    extinction_threshold = 20 # FIXME
    pop_size = 3000000 # FIXME

    if len(display_scenario_names) < len(scenario_names):
        display_scenario_names = scenario_names

    all_last_days_with_infections_exclude_extinction = []

    all_herd_immunity_thresholds_exclude_extinction = []
    all_herd_immunity_thresholds_day_exclude_extinction = []

    for scenario in scenario_names:
        with multiprocessing.Pool(processes=4) as pool:
            # Get mean Rt per day
            mean_rt_by_day_exclude_extinction = pool.starmap(get_effective_r_by_day,
                                                [(output_dir, scenario, exp_id, num_days, extinction_threshold) for exp_id in experiment_ids])
            mean_rt_by_day_exclude_extinction = [x for x in mean_rt_by_day_exclude_extinction if x is not None] # Remove extinction cases
            # Get herd immunity threshold = proportion of population no longer susceptible on last day for which Rt >= 1
            herd_immunity_threshold_exclude_extinction = pool.starmap(get_herd_immunity_threshold,
                                                [(output_dir, scenario, exp_id, num_days, pop_size, extinction_threshold) for exp_id in experiment_ids])
            herd_immunity_threshold_exclude_extinction = [run for run in herd_immunity_threshold_exclude_extinction if not np.isnan(run)]
            # Get day on which herd immunity threshold is reached
            herd_immunity_threshold_day_exclude_extinction = pool.starmap(get_herd_immunity_threshold,
                                                [(output_dir, scenario, exp_id, num_days, pop_size, extinction_threshold, True) for exp_id in experiment_ids])
            herd_immunity_threshold_day_exclude_extinction = [run for run in herd_immunity_threshold_day_exclude_extinction if not np.isnan(run)]

            all_herd_immunity_thresholds_exclude_extinction.append(herd_immunity_threshold_exclude_extinction)
            all_herd_immunity_thresholds_day_exclude_extinction.append(herd_immunity_threshold_day_exclude_extinction)

            # Calculate last day with new infections
            last_day_with_infections = []
            for run in new_cases_per_day:
                for day in range(num_days - 1, -1, -1):
                    if run[day] > 0:
                        last_day_with_infections.append(day)
                    elif day == 0:
                        last_day_with_infections.append(day)
            last_day_with_infections_exclude_extinction = [x for x in last_day_with_infections if x > 35]
            all_last_days_with_infections_exclude_extinction.append(last_day_with_infections_exclude_extinction)

            # Plot evolution of new cases per day
            plot_evolution(output_dir, "evolution", scenario, new_cases_per_day, num_days, 170000)

            # Plot evolution of Rt
            plot_effective_r_by_day(output_dir, "mean_rt_exclude_extinction_avg", scenario, mean_rt_by_day_exclude_extinction, num_days)


    plot_extinction_probabilities(output_dir, "extinction_probabilities", display_scenario_names, all_total_cases, extinction_threshold)
    plot_ar(output_dir, "ar_exclude_extinction", display_scenario_names, all_total_cases, num_days, pop_size, extinction_threshold)
    plot_ar(output_dir, "ar_exclude_extinction", display_scenario_names, all_total_cases, num_days, pop_size, 0, True)
    plot_day_last_infection(output_dir, "day_of_last_infection_exclude_extinction", display_scenario_names, all_last_days_with_infections_exclude_extinction)

    plot_herd_immunity_threshold(output_dir, "herd_immunity_threshold_exclude_extinction_smoothed_rt", display_scenario_names, all_herd_immunity_thresholds_exclude_extinction)
    plot_herd_immunity_threshold(output_dir, "herd_immunity_threshold_day_exclude_extinction_smoothed_rt", display_scenario_names, all_herd_immunity_thresholds_day_exclude_extinction, show_day=True)
'''

import argparse
import matplotlib.pyplot as plt
import multiprocessing

from plots import plot_p80s, plot_final_size_frequencies, plot_new_cases_per_day
from util import get_experiment_ids, get_cases_output, get_p80, get_total_cases

def main(output_dir, scenario_names, display_scenario_names):
    num_days = 200

    all_p80s = []
    all_final_sizes = []
    for s_i in range(len(scenario_names)):
        scenario_name = scenario_names[s_i]
        print(scenario_name)

        experiment_ids = get_experiment_ids(output_dir, scenario_name)
        with multiprocessing.Pool(processes=4) as pool:
            cases_output = pool.starmap(get_cases_output, [(output_dir, scenario_name, exp_id) for exp_id in experiment_ids])
            # Calculate the proportion of infected individuals responsible for 80% of infections
            all_p80s.append([get_p80(output["secondary_cases_by_individual"]) for output in cases_output]) # TODO exclude extinction?
            # Get final outbreak size
            total_cases = [get_total_cases(output["cases_by_day"], num_days) for output in cases_output]
            all_final_sizes.append(total_cases)

            # Sort total cases from high to low & print
            # Used to determine extinction threshold
            total_cases.sort(reverse=True)
            print(total_cases)

            plot_new_cases_per_day(output_dir, "new_cases_per_day", display_scenario_names, [output["cases_by_day"] for output in cases_output], num_days)
    plot_p80s(output_dir, "p80s", display_scenario_names, all_p80s)
    plot_final_size_frequencies(output_dir, "final_size_frequencies", display_scenario_names, all_final_sizes, num_days)

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", type=str, help="Directory containing simulation results")
    parser.add_argument("scenario_names", type=str, nargs="+", help="Names of scenarios to be postprocessed")
    parser.add_argument("--display_scenario_names", type=str, nargs="+", default=[], help="Names for scenarios to be displayed on plots")

    args = parser.parse_args()
    main(args.output_dir, args.scenario_names, args.display_scenario_names)
