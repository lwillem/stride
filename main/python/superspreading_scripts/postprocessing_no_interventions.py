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

import argparse
import matplotlib.pyplot as plt
import multiprocessing

from plots import plot_p80s, plot_final_size_frequencies, plot_day_last_infection, plot_effective_r_by_day
from plots import plot_ar, plot_cumulative_cases_per_day, plot_new_cases_per_day, plot_extinction_probabilities, plot_herd_immunity_threshold

from util import get_experiment_ids, get_cases_output, get_p80, get_total_cases, get_day_of_last_infection, get_rt_by_day, get_herd_immunity_threshold

def main(output_dir, scenario_names, display_scenario_names):

    if len(display_scenario_names) < len(scenario_names):
        display_scenario_names = scenario_names

    num_days = 200
    extinction_threshold = 20
    population_size = 3000000

    all_days_last_infection = []
    all_p80s = []
    all_final_sizes = []
    all_hits = []
    all_hits_day = []

    for s_i in range(len(scenario_names)):
        scenario_name = scenario_names[s_i]
        print(scenario_name)

        experiment_ids = get_experiment_ids(output_dir, scenario_name)
        with multiprocessing.Pool(processes=4) as pool:
            cases_output = pool.starmap(get_cases_output, [(output_dir, scenario_name, exp_id) for exp_id in experiment_ids])
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
            all_hits_day.append([get_herd_immunity_threshold(output["rt_by_day"], output["cases_by_day"], num_days, output["parameters"]["population_size"], get_day=True) for output in cases_output])

            # Sort total cases from high to low & print
            # Used to determine extinction threshold
            total_cases.sort(reverse=True)
            print(total_cases)

            plot_new_cases_per_day(output_dir, "new_cases_per_day", scenario_name, [output["cases_by_day"] for output in cases_output], num_days)
            plot_cumulative_cases_per_day(output_dir, "cumulative_cases_per_day", scenario_name, [output["cases_by_day"] for output in cases_output], num_days)
            plot_effective_r_by_day(output_dir, "rt_by_day", scenario_name, [output["rt_by_day"] for output in cases_output], num_days)

    plot_extinction_probabilities(output_dir, "extinction_probabilities", display_scenario_names, all_final_sizes, extinction_threshold)
    plot_p80s(output_dir, "p80s", display_scenario_names, all_p80s)
    plot_ar(output_dir, "ar", display_scenario_names, all_final_sizes, num_days, population_size)
    plot_ar(output_dir, "ar_exclude_extinction", display_scenario_names, all_final_sizes, num_days, population_size, extinction_threshold=extinction_threshold)
    plot_herd_immunity_threshold(output_dir, "hits", display_scenario_names, all_hits)
    plot_herd_immunity_threshold(output_dir, "hits_day", display_scenario_names, all_hits_day, show_day=True)
    plot_final_size_frequencies(output_dir, "final_size_frequencies", display_scenario_names, all_final_sizes, num_days)
    plot_day_last_infection(output_dir, "last_day_with_infections", display_scenario_names, all_days_last_infection, num_days)

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", type=str, help="Directory containing simulation results")
    parser.add_argument("scenario_names", type=str, nargs="+", help="Names of scenarios to be postprocessed")
    parser.add_argument("--display_scenario_names", type=str, nargs="+", default=[], help="Names for scenarios to be displayed on plots")

    args = parser.parse_args()
    main(args.output_dir, args.scenario_names, args.display_scenario_names)
