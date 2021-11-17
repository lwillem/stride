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


    """

    """
    """

    import argparse
    import multiprocessing

    from plots import plot_ar, plot_cumulative_cases_per_day, plot_day_of_last_infection, plot_day_of_peak, plot_extinction_probabilities, plot_effective_r_by_day, plot_final_size_frequencies, plot_herd_immunity_threshold, plot_new_cases_per_day, plot_p80s, plot_peak_sizes, plot_secondary_cases_distribution, plot_transmissions_by_location

    from util import get_day_of_last_infection, get_experiment_ids, get_herd_immunity_threshold, get_output, get_p80, get_total_cases

    def main(output_dir):
        baseline_scenario_name = "baseline"

        overdispersion_scenario_names = ["infectiousness_overdispersion", "contacts_overdispersion"]
        overdispersion_parameters = ["1000", "100", "60", "40", "20"]

        num_days = 200
        population_size = 3000000

        extinction_threshold = 20

        for overdispersion_scenario in overdispersion_scenario_names:
            scenario_names = [baseline_scenario_name] + [overdispersion_scenario + "_" + overdispersion for overdispersion in overdispersion_parameters]

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

            all_secondary_cases = []

            alpha = r"$\alpha$"
            display_scenario_names = ["Baseline", alpha + " = 10", alpha + " = 1", alpha + " = 0.6", alpha + " = 0.4", alpha + " = 0.2"]

            for scenario_name in scenario_names:
                print(scenario_name)
                experiment_ids = get_experiment_ids(output_dir, scenario_name)
                output = {}
                with multiprocessing.Pool(processes=4) as pool:
                    output = pool.starmap(get_output, [(output_dir, scenario_name, exp_id) for exp_id in experiment_ids])

                    # Sort total cases from high to low & print
                    # Used to determine extinction threshold
                    total_cases = [get_total_cases(run["cases_per_day"], num_days) for run in output]
                    total_cases.sort(reverse=True)
                    print(total_cases)

                    all_final_sizes.append(total_cases)
                    all_cases_per_day.append([run["cases_per_day"] for run in output])

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

                    plot_new_cases_per_day(output_dir, "new_cases_per_day", scenario_name, [run["cases_per_day"] for run in output], num_days, y_max=170000)
                    plot_cumulative_cases_per_day(output_dir, "cumulative_cases_per_day", scenario_name, [run["cases_per_day"] for run in output], num_days, y_max=3000000)

                    plot_effective_r_by_day(output_dir, "rt", scenario_name, [run["rt_by_day"] for run in output], num_days, y_max=30)

                    plot_transmissions_by_location(output_dir, "transmissions_by_location", scenario_name, [run["transmissions_by_location"] for run in output])

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

            plot_secondary_cases_distribution(output_dir, "secondary_cases_distribution_" + overdispersion_scenario, display_scenario_names, all_secondary_cases)


    if __name__=="__main__":
        parser = argparse.ArgumentParser()
        parser.add_argument("output_dir", type=str, help="Directory containing simulation results")

        args = parser.parse_args()
        main(args.output_dir)

    """
