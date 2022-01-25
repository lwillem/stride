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
    Script to create plots comparing theoretical description estimates
    to simulation results.
"""
"""
import csv
import matplotlib.pyplot as plt
import numpy as np
import os

from estimate_transmission_probability import estimate_effective_contacts
from util import save_figure

def main(output_dir, scenario_names, overdispersion_params, population_file, contact_matrix_file):
    infectious_period_length = 7 # FIXME This probably should not be hard-coded here...

    for scenario_i in range(len(scenario_names)):

        # Estimate mean + variance number of effective contacts for index cases

        means_theoretical_by_tp, var_theoretical_by_tp = estimate_effective_contacts(population_file, contact_matrix_file,
                                                                tps_sorted, infectious_period_length,
                                                                overdispersion=overdispersion_vals[scenario_i],
                                                                person_ids=person_ids)



if __name__=="__main__":
    parser.add_argument("--contact_matrix_file", type=str, default=os.path.join("..", "resources", "data", "contact_matrix_flanders_conditional_teachers.xml"))
"""

import argparse
import multiprocessing

from plots import plot_comparison_means_theoretical_sims, plot_comparison_variance_theoretical_sims

from postprocessing_util import get_experiment_ids, get_num_secondary_cases_per_index_case

def main(output_dir, index_case_id, num_parallel_workers):
    baseline_scenario_name = "simplified_baseline"
    overdispersion_scenario_names = ["infectiousness", "contacts"]
    overdispersion_parameters = ["1000", "100", "60", "40", "20"]
    transmission_probabilities = [0.025, 0.05, 0.075, 0.1]

    for overdispersion_scenario in overdispersion_scenario_names:
        scenario_names = [baseline_scenario_name] + ["simplified_" + overdispersion_scenario + "_overdispersion_" + overdispersion for overdispersion in overdispersion_parameters]
        for scenario_name in scenario_names:
            print(scenario_name)

            secondary_cases_by_tp = []

            for tp in transmission_probabilities:
                full_scenario_name = scenario_name + "_pid_" + str(index_case_id) + "_tp_" + str(tp)
                experiment_ids = get_experiment_ids(output_dir, full_scenario_name)
                with multiprocessing.Pool(processes=num_parallel_workers) as pool:
                    secondary_cases = pool.starmap(get_num_secondary_cases_per_index_case, [(output_dir, full_scenario_name, exp_id) for exp_id in experiment_ids])
                    secondary_cases_by_tp.append(secondary_cases)

            theoretical_means = secondary_cases_by_tp # TODO
            theoretical_variances = [0] * len(transmission_probabilities) # TODO

            plot_comparison_means_theoretical_sims(output_dir, "mean_comparison_", scenario_name + "_pid_" + str(index_case_id), transmission_probabilities, secondary_cases_by_tp, theoretical_means)
            plot_comparison_variance_theoretical_sims(output_dir, "variance_comparison_", scenario_name + "_pid_" + str(index_case_id), transmission_probabilities, secondary_cases_by_tp, theoretical_variances)


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", type=str, help="Directory containing simulation results")
    parser.add_argument("index_case_id", type=int, help="Id of index case.")
    parser.add_argument("--num_parallel_workers", type=int, default=4)

    args = parser.parse_args()
    main(args.output_dir, args.index_case_id, args.num_parallel_workers)
