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
#  Copyright 2021, Kuylen E
############################################################################ #

"""
"""

import argparse
import multiprocessing
import numpy as np

from plots import plot_resurgence_probabilities
from plots import plot_ar, plot_cumulative_cases_per_day, plot_day_of_last_infection, plot_extinction_probabilities, plot_new_cases_per_day, plot_num_cases_over_period, plot_effective_r_by_day, plot_transmissions_by_location, plot_final_size_frequencies, plot_peak_sizes, plot_p80s, plot_secondary_cases_distribution, plot_herd_immunity_threshold
from postprocessing_util import get_num_secondary_cases_frequencies, get_herd_immunity_threshold, get_cases_over_period, get_experiment_ids, get_output_per_day, get_summary_output, get_total_cases, get_transmissions_by_location, get_day_of_last_infection, get_resurgence_probability

def main(output_dir, num_parallel_workers):
    num_days = np.nan
    population_size = np.nan

    extinction_threshold = 10000
    resurgence_threshold = 500

    baseline_scenario_name = "sd_baseline"

    overdispersion_scenario_names = ["sd_infectiousness_overdispersion", "sd_contacts_overdispersion"]
    overdispersion_parameters = ["1000", "100", "60", "40", "20"]

    alpha = r"$\alpha$"
    display_scenario_names = ["Baseline", alpha + " = 10", alpha + " = 1", alpha + " = 0.6", alpha + " = 0.4", alpha + " = 0.2"]

    events = {
        "Start lockdown": 30,
        "End lockdown": 90
    }

    num_days = np.nan
    population_size = np.nan

    for overdispersion_scenario in overdispersion_scenario_names:
        scenario_names = [baseline_scenario_name] + [overdispersion_scenario + "_" + overdispersion for overdispersion in overdispersion_parameters]

        all_cases_per_day = []
        all_cases_before_lockdown = []
        all_cases_during_lockdown = []
        all_cases_after_lockdown = []

        all_final_sizes = []
        all_hits = []

        all_p80s = []
        all_secondary_cases = []

        resurgence_probabilities = []
        for scenario_name in scenario_names:
            print(scenario_name)
            experiment_ids = get_experiment_ids(output_dir, scenario_name)

            with multiprocessing.Pool(processes=num_parallel_workers) as pool:
                output_per_day = pool.starmap(get_output_per_day, [(output_dir, scenario_name, exp_id) for exp_id in experiment_ids])
                output_summary = pool.starmap(get_summary_output, [(output_dir, scenario_name, exp_id) for exp_id in experiment_ids])
                transmissions_by_location = pool.starmap(get_transmissions_by_location, [(output_dir, scenario_name, exp_id) for exp_id in experiment_ids])
                output_num_secondary_cases = pool.starmap(get_num_secondary_cases_frequencies, [(output_dir, scenario_name, exp_id) for exp_id in experiment_ids])

                num_days = output_summary[0]["num_days"]
                population_size = output_summary[0]["population_size"]

                resurgence_probabilities.append(get_resurgence_probability([run["cases_per_day"] for run in output_per_day], 30, 90, 600, 500))

                cases_before_lockdown = [get_cases_over_period(run["cases_per_day"], 0, 30) for run in output_per_day]
                cases_during_lockdown = [get_cases_over_period(run["cases_per_day"], 30, 90) for run in output_per_day]
                cases_after_lockdown = [get_cases_over_period(run["cases_per_day"], 90, 600) for run in output_per_day]

                #cases_before_lockdown.sort(reverse=True)
                #print("CASES BEFORE LOCKDOWN")
                #print(cases_before_lockdown)

                #cases_during_lockdown.sort(reverse=True)
                #print("CASES DURING LOCKDOWN")
                #print(cases_during_lockdown)

                #cases_after_lockdown.sort(reverse=True)
                #print("CASES RELEASE")
                #print(cases_after_lockdown)

                all_cases_per_day.append([run["cases_per_day"] for run in output_per_day])

                total_cases = [get_total_cases(run["cases_per_day"], num_days) for run in output_per_day]
                #total_cases.sort(reverse=True)
                #print("TOTAL CASES")
                #print(total_cases)

                all_cases_before_lockdown.append(cases_before_lockdown)
                all_cases_during_lockdown.append(cases_during_lockdown)
                all_cases_after_lockdown.append(cases_after_lockdown)

                all_final_sizes.append(total_cases)
                all_hits.append([get_herd_immunity_threshold(run["rt_by_day"], run["cases_per_day"], num_days, population_size) for run in output_per_day])

                all_p80s.append([run["p80"] for run in output_summary])
                all_secondary_cases.append(output_num_secondary_cases)

                plot_cumulative_cases_per_day(output_dir, "cumulative_cases_per_day", scenario_name, [run["cases_per_day"] for run in output_per_day], num_days, y_max=population_size, events=events)
                plot_new_cases_per_day(output_dir, "new_cases_per_day", scenario_name, [run["cases_per_day"] for run in output_per_day], num_days, y_max=250000, events=events)

                plot_effective_r_by_day(output_dir, "rt", scenario_name, [run["rt_by_day"] for run in output_per_day], num_days, y_max=20, events=events)
                plot_transmissions_by_location(output_dir, "transmissions_by_location", scenario_name, transmissions_by_location)

        plot_ar(output_dir, "ar_" + overdispersion_scenario, display_scenario_names, all_final_sizes, num_days, population_size, violin_plot=True)

        plot_num_cases_over_period(output_dir, "num_cases_before_lockdown_" + overdispersion_scenario, display_scenario_names, 0, 30, all_cases_before_lockdown,y_min=-100, y_max=55000)
        plot_num_cases_over_period(output_dir, "num_cases_during_lockdown_" + overdispersion_scenario, display_scenario_names, 30, 90, all_cases_during_lockdown, y_min=-100, y_max=95000)
        plot_num_cases_over_period(output_dir, "num_cases_after_lockdown_" + overdispersion_scenario, display_scenario_names, 90, 600, all_cases_after_lockdown,y_min=-20000, y_max=7000000)
        plot_num_cases_over_period(output_dir, "num_cases_after_lockdown_exclude_extinction" + overdispersion_scenario, display_scenario_names, 90, 600, all_cases_after_lockdown,y_min=-20000, y_max=7000000, extinction_threshold=resurgence_threshold)

        plot_herd_immunity_threshold(output_dir, "hit_" + overdispersion_scenario, display_scenario_names, all_hits, show_day=False)
        plot_peak_sizes(output_dir, "peak_sizes_" + overdispersion_scenario, display_scenario_names, all_cases_per_day, 90, num_days, ymin=-500, ymax=300000)
        plot_peak_sizes(output_dir, "peak_sizes_exclude_no_resurgence_" + overdispersion_scenario, display_scenario_names, all_cases_per_day, 90, num_days, ymin=-500, ymax=300000, extinction_threshold=500)

        plot_p80s(output_dir, "p80s_" + overdispersion_scenario, display_scenario_names, all_p80s)
        plot_secondary_cases_distribution(output_dir, "secondary_cases_distribution_" + overdispersion_scenario, display_scenario_names, all_secondary_cases)

        plot_resurgence_probabilities(output_dir, "resurgence_probabilities_" + overdispersion_scenario, display_scenario_names, resurgence_probabilities)


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", type=str)
    parser.add_argument("--num_parallel_workers", type=int, default=4)

    args = parser.parse_args()
    main(args.output_dir, args.num_parallel_workers)
