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
import argparse
import csv
import matplotlib.pyplot as plt
import multiprocessing
import numpy as np
import os

from estimate_transmission_probability import estimate_effective_contacts
from util import get_trans_prob_by_exp, save_figure

def get_index_case_ids(output_dir, scenario_name, experiment_id):
    index_case_ids = []

    transmissions_file = os.path.join(output_dir, scenario_name, "exp" + "{:04}".format(experiment_id), "event_log.txt")
    with open(transmissions_file) as f:
        for line in f:
            line = line.split(" ")
            tag = line[0]
            if tag == "[PRIM]":
                index_case_id = int(float(line[1]))
                index_case_ids.append(index_case_id)

    return (experiment_id, index_case_ids)

def get_secondary_cases_per_index_case(output_dir, scenario_name, experiment_id):
    secondary_cases = {}

    transmissions_file = os.path.join(output_dir, scenario_name, "exp" + "{:04}".format(experiment_id), "event_log.txt")
    with open(transmissions_file) as f:
        for line in f:
            line = line.split(" ")
            tag = line[0]
            if tag == "[PRIM]":
                index_case_id = int(float(line[1]))
                secondary_cases[index_case_id] = 0
            elif tag == "[TRAN]":
                infector_id = int(float(line[2]))
                if infector_id in secondary_cases:
                    secondary_cases[infector_id] += 1
                else:
                    secondary_cases[infector_id] = 1
    if len(secondary_cases) > 1:
        print("WARNING: more than 1 index case")
    secondary_cases_per_index_case = sum(secondary_cases.values()) / len(secondary_cases)

    return (experiment_id, secondary_cases_per_index_case)

def main(output_dir, scenario_names, overdispersion_params, population_file, contact_matrix_file):
    infectious_period_length = 7 # FIXME This probably should not be hard-coded here...

    for scenario_i in range(len(scenario_names)):
        print(scenario_names[scenario_i])

        # Get transmission probability per experiment
        experiments = get_trans_prob_by_exp(output_dir, scenario)

        secondary_cases = []
        index_case_ids = []
        with multiprocessing.Pool(processes=4) as pool:
            secondary_cases = pool.starmap(get_secondary_cases_per_index_case,
                                            [(output_dir, scenario_names[scenario_i], exp_id) for exp_id in experiments.keys()])
            index_case_ids = pool.starmap(get_index_case_ids,
                                            [(output_dir, scenario_names[scenario_i], exp_id) for exp_id in experiments.keys()])

        # Group by transmission probability
        secondary_cases_by_tp = {}
        for experiment_id, cases in secondary_cases:
            tp = experiments[experiment_id]
            if tp in secondary_cases_by_tp:
                secondary_cases_by_tp[tp].append(cases)
            else:
                secondary_cases_by_tp[tp] = [cases]

        tps_sorted = list(secondary_cases_by_tp.keys())
        tps_sorted.sort()

        # Estimate mean + variance number of effective contacts for index cases
        person_ids = list(set([x[1][0] for x in index_case_ids]))

        means_theoretical_by_tp, var_theoretical_by_tp = estimate_effective_contacts(population_file, contact_matrix_file,
                                                                tps_sorted, infectious_period_length,
                                                                overdispersion=overdispersion_vals[scenario_i],
                                                                person_ids=person_ids)

        # Plot theoretical estimate of mean number of secondary cases per index case
        # VS mean and 95% interval of number of secondary cases per index case from simulations
        means = [np.mean(secondary_cases_by_tp[tp]) for tp in tps_sorted]
        lower = [np.percentile(secondary_cases_by_tp[tp], 2.5) for tp in tps_sorted]
        upper = [np.percentile(secondary_cases_by_tp[tp], 97.5) for tp in tps_sorted]

        plt.plot(tps_sorted, means, marker="o", label="Simulations")
        plt.fill_between(tps_sorted, lower, upper, color="lightgrey")

        plt.plot(tps_sorted, [means_theoretical_by_tp[tp] for tp in tps_sorted], linestyle="None", marker="^", label="Theoretical")

        plt.xlabel("Mean transmission probality")
        plt.ylabel("Secondary cases caused by index case")
        plt.ylim(-2, 100)
        plt.legend()
        save_figure(output_dir, "mean_ibrn_comparison_" + scenario_names[scenario_i], extension="png")

        # Plot theoretical estimate of variance for number of secondary cases per index case
        # VS variance for number of secondary cases per index case from simulations

        variances = [np.var(secondary_cases_by_tp[tp]) for tp in tps_sorted]

        plt.plot(tps_sorted, variances, linestyle="None", marker="o")
        plt.plot(tps_sorted, [var_theoretical_by_tp[tp] for tp in tps_sorted], linestyle="None", marker="^")
        plt.legend(["Simulations", "Theoretical"])
        plt.xlabel("Mean transmission probality")
        plt.ylabel("Variance of secondary cases caused by index case")
        save_figure(output_dir, "var_" + scenario_names[scenario_i], extension="png")


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", type=str, help="Directory containing simulation results")
    parser.add_argument("scenario_names", type=str, nargs="+", help="Names of scenarios to be postprocessed")
    parser.add_argument("--overdispersion_params", type=float, nargs="+", default=[None, 10, 0.6, 0.4, 0.2], help="Overdispersion parameters for the scenarios - Use None if not applicable")
    parser.add_argument("--population_file", type=str, default=os.path.join("..", "resources", "data", "pop_belgium3000k_c500_teachers_censushh.csv"))
    parser.add_argument("--contact_matrix_file", type=str, default=os.path.join("..", "resources", "data", "contact_matrix_flanders_conditional_teachers.xml"))

    args = parser.parse_args()
    main(args.output_dir, args.scenario_names, args.overdispersion_params, args.population_file, args.contact_matrix_file)"""
