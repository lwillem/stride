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
    Postprocessing to get mean number of secondary cases per index case
    + comparison of offspring distribution to negative binomial distribution.
    Using results from simulations tracking only secondary cases from index case.
"""
"""
import argparse
import multiprocessing
import matplotlib.pyplot as plt
import numpy as np
import os

from scipy.stats import nbinom, probplot

from util import get_experiment_ids, get_trans_prob_by_exp, save_figure

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
                    print("Not an index case")
    if len(secondary_cases) > 1:
        print("WARNING: more than 1 index case")

    secondary_cases_per_index_case = sum(secondary_cases.values()) / len(secondary_cases)
    return (experiment_id, secondary_cases_per_index_case)

def main(output_dir, scenario_names, overdispersion_params, display_scenario_names):
    if len(display_scenario_names) < len(scenario_names):
        display_scenario_names = scenario_names

    all_means = []
    all_means_exclude_extinction = []

    for scenario_i in range(len(scenario_names)):
        print(scenario_names[scenario_i])

        experiments = get_trans_prob_by_exp(output_dir, scenario)

        secondary_cases = []
        with multiprocessing.Pool(processes=4) as pool:
            secondary_cases = pool.starmap(get_secondary_cases_per_index_case,
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

        # Mean secondary cases per index case
        mean_secondary_cases_by_tp = [np.mean(secondary_cases_by_tp[tp]) for tp in tps_sorted]
        all_means.append(mean_secondary_cases_by_tp)

        # Mean secondary cases per index case,
        # excluding runs where index case makes 0 secondary cases
        mean_secondary_cases_by_tp_exclude_extinction = []
        for tp in tps_sorted:
            if tp == 0:
                mean_secondary_cases_by_tp_exclude_extinction.append(np.nan)
            else:
                mean_secondary_cases_by_tp_exclude_extinction.append(np.mean([x for x in secondary_cases_by_tp[tp] if x > 0]))
        all_means_exclude_extinction.append(mean_secondary_cases_by_tp_exclude_extinction)

        # Create QQ-plots to compare offspring distribution
        # to negative binomial distribution
        tp_i = 0
        for tp in tps_sorted:
            k = overdispersion_params[scenario_i]
            if k is None:
                k = np.inf
            res = probplot(secondary_cases_by_tp[tp], dist=nbinom, sparams=(k, k / (mean_secondary_cases_by_tp[tp_i] + k)), fit=False, plot=plt)
            save_figure(output_dir, "QQplot_" + scenario + "_tp_" + str(tp))

            tp_i += 1

    # Plot mean secondary cases per index case
    for scenario_i in scenario_names:
        plt.plot([0.0, 0.025, 0.05, 0.075, 0.10], all_means[scenario_i])
    plt.xlabel("Mean individual transmission probability")
    plt.xticks([0.00, 0.02, 0.04, 0.06, 0.08, 0.10])
    plt.ylabel("Mean number of secondary cases per index case")
    plt.legend(display_scenario_names)
    save_figure(output_dir, "mean_secondary_cases_by_tp")

    # Plot mean secondary cases per index case, excluding runs where index case had 0 secondary cases
    for scenario_i in scenario_names:
        plt.plot([0.0, 0.025, 0.05, 0.075, 0.10], all_means_exclude_extinction[scenario_i])
    plt.xlabel("Mean individual transmission probability")
    plt.xticks([0.00, 0.02, 0.04, 0.06, 0.08, 0.10])
    plt.ylabel("Mean number of secondary cases per index case")
    plt.legend(display_scenario_names)
    save_figure(output_dir, "mean_secondary_cases_by_tp_exclude_extinction")


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", type=str, help="Directory containing simulation results")
    parser.add_argument("scenario_names", type=str, nargs="+", help="Names of scenarios to be postprocessed")
    parser.add_argument("--overdispersion_params", type=float, nargs="+", default=[None, 10,0.6,0.4,0.2,] help="Overdispersion parameters for the scenarios")
    parser.add_argument("--display_scenario_names", type=str, nargs="+", default=[], help="Names for scenarios to be displayed on plots")

    args = parser.parse_args()
    main(args.output_dir, args.scenario_names, args.overdispersion_params, args.display_scenario_names)"""
