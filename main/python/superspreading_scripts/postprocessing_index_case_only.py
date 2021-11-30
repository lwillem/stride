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
def main(output_dir):

    for overdispersion_scenario in overdispersion_scenario_names:

        for scenario_name in scenario_names:

            with multiprocessing.Pool(processes=4) as pool:



"""

import argparse
import multiprocessing

from plots import plot_qq, plot_offspring_distributions, plot_secondary_cases_per_index_case
from postprocessing_util import get_experiment_ids, get_num_secondary_cases_per_index_case, get_summary_output

def main(output_dir, num_parallel_workers):
    baseline_scenario_name = "ico_baseline"
    overdispersion_scenario_names = ["ico_infectiousness_overdispersion", "ico_contacts_overdispersion"]
    overdispersion_parameters = ["1000", "100", "60", "40", "20"]

    transmission_probabilities = [0.025, 0.05, 0.075, 0.1]

    alpha = r"$\alpha$"
    display_scenario_names = ["Baseline", alpha + " = 10", alpha + " = 1", alpha + " = 0.6", alpha + " = 0.4", alpha + " = 0.2"]

    for overdispersion_scenario in overdispersion_scenario_names:
        all_secondary_cases_by_tp = []

        scenario_names = [baseline_scenario_name] + [overdispersion_scenario + "_" + overdispersion for overdispersion in overdispersion_parameters]

        for s_i in range(len(scenario_names)):
            scenario_name = scenario_names[s_i]
            print(scenario_name)

            secondary_cases_by_tp = {}

            for tp in transmission_probabilities:
                full_scenario_name = scenario_name + "_tp_" + str(tp)
                experiment_ids = get_experiment_ids(output_dir, full_scenario_name)

                with multiprocessing.Pool(processes=num_parallel_workers) as pool:
                    summary_output = pool.starmap(get_summary_output, [(output_dir, full_scenario_name, exp_id) for exp_id in experiment_ids])
                    secondary_cases_per_index_case = pool.starmap(get_num_secondary_cases_per_index_case, [(output_dir, full_scenario_name, exp_id) for exp_id in experiment_ids])

                    secondary_cases_by_tp[tp] = secondary_cases_per_index_case

            all_secondary_cases_by_tp.append(secondary_cases_by_tp)

            plot_offspring_distributions(output_dir, "offspring_distribution", scenario_name, secondary_cases_by_tp)
            if s_i != 0:
                k = int(overdispersion_parameters[s_i - 1]) / 100
                plot_qq(output_dir, "qq_plot", scenario_name, k, secondary_cases_by_tp)

        plot_secondary_cases_per_index_case(output_dir, "secondary_cases_per_index_case_" + overdispersion_scenario, display_scenario_names, all_secondary_cases_by_tp)
        plot_secondary_cases_per_index_case(output_dir, "secondary_cases_per_index_case_exclude_extinction" + overdispersion_scenario, display_scenario_names, all_secondary_cases_by_tp, exclude_extinction=True)

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", type=str, help="Directory containing simulation results")
    parser.add_argument("--num_parallel_workers", type=int, default=4)

    args = parser.parse_args()
    main(args.output_dir, args.num_parallel_workers)
