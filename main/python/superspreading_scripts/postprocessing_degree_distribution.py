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
import matplotlib.pyplot as plt
import multiprocessing

from plots import save_figure
from postprocessing_util import get_degree_distribution_from_file, get_experiment_ids, get_summary_output

def main(output_dir, num_parallel_workers):
    scenario_names = ["dd_baseline",
                        "dd_contacts_overdispersion_1000",
                        "dd_contacts_overdispersion_100",
                        "dd_contacts_overdispersion_60",
                        "dd_contacts_overdispersion_40",
                        "dd_contacts_overdispersion_20"]
    alpha = r"$\alpha$"
    display_scenario_names = ["Baseline", alpha + " = 10", alpha + " = 1", alpha + " = 0.6", alpha + " = 0.4", alpha + " = 0.2"]

    for scenario_name in scenario_names:
        scenario_name = output_prefix + scenario_name
        print(scenario_name)
        experiment_ids = get_experiment_ids(output_dir, scenario_name)

        with multiprocessing.Pool(processes=num_parallel_workers) as pool:
            degree_distribution = pool.starmap(get_degree_distribution_from_file, [(output_dir, scenario_name, exp_id) for exp_id in experiment_ids])
            summary_output = pool.starmap(get_summary_output, [(output_dir, scenario_name, exp_id, "summary.csv") for exp_id in experiment_ids])
            population_size = summary_output[0]["population_size"]

            for run in degree_distribution:
                degrees = list(run.keys())
                degrees.sort()
                plt.plot(degrees, [run[d] / population_size for d in degrees])
    plt.legend(display_scenario_names)
    plt.xlim(0, 250)
    plt.xlabel("Degree")
    plt.ylabel("Frequency")

    save_figure(output_dir, "degree_distributions")

if __name__=="__main__":
     parser = argparse.ArgumentParser()

     parser.add_argument("output_dir", type=str, help="Directory containing simulation results")
     parser.add_argument("output_prefix", type=str, default="")
     parser.add_argument("--num_parallel_workers", type=int, default=4)

     args = parser.parse_args()

     main(args.output_dir, args.num_parallel_workers)
