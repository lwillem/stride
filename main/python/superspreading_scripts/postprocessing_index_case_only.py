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

import argparse
import multiprocessing

from plots import plot_qq, plot_offspring_distributions, plot_secondary_cases_per_index_case
from util import get_experiment_ids, get_output

def main(output_dir):
    baseline_scenario_name = "ico_baseline"
    overdispersion_scenario_names = ["ico_infectiousness_overdispersion", "ico_contacts_overdispersion"]
    overdispersion_parameters = ["1000", "100", "60", "40", "20"]

    num_days = 40
    population_size = 3000000

    for overdispersion_scenario in overdispersion_scenario_names:
        scenario_names = [baseline_scenario_name] + [overdispersion_scenario + "_" + overdispersion for overdispersion in overdispersion_parameters]

        alpha = r"$\alpha$"
        display_scenario_names = ["Baseline", alpha + " = 10", alpha + " = 1", alpha + " = 0.6", alpha + " = 0.4", alpha + " = 0.2"]

        all_secondary_cases_by_tp = []
        i = 0
        for scenario_name in scenario_names:
            print(scenario_name)
            experiment_ids = get_experiment_ids(output_dir, scenario_name)

            output = {}
            with multiprocessing.Pool(processes=4) as pool:
                output = pool.starmap(get_output, [(output_dir, scenario_name, exp_id) for exp_id in experiment_ids])

                # Group number of secondary cases / individual by transmission probability
                secondary_cases_by_tp = {}
                for run in output:
                    tp = run["parameters"]["transmission_probability"]
                    if tp in secondary_cases_by_tp:
                        secondary_cases_by_tp[tp] += list(run["secondary_cases_by_index_case"].values())
                    else:
                        secondary_cases_by_tp[tp] = list(run["secondary_cases_by_index_case"].values())

                all_secondary_cases_by_tp.append(secondary_cases_by_tp)

                if i != 0:
                    k = int(overdispersion_parameters[i - 1]) / 100
                    plot_qq(output_dir, "qq_plot", scenario_name, k, secondary_cases_by_tp)

                plot_offspring_distributions(output_dir, "offspring_distribution", scenario_name, secondary_cases_by_tp)


            i += 1

        plot_secondary_cases_per_index_case(output_dir, "secondary_cases_per_index_case_" + overdispersion_scenario, display_scenario_names, all_secondary_cases_by_tp)
        plot_secondary_cases_per_index_case(output_dir, "secondary_cases_per_index_case_exclude_extinction" + overdispersion_scenario, display_scenario_names, all_secondary_cases_by_tp, exclude_extinction=True)


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", type=str, help="Directory containing simulation results")

    args = parser.parse_args()
    main(args.output_dir)
