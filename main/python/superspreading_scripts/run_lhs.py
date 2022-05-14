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
#  Copyright 2022, Kuylen E
############################################################################ #

import argparse

from scipy.stats import qmc

from run_util import get_mean_non_truncated_gamma, run_parallel

def run_lhs(disease_config_file, holidays_file, mean_transmission_probability, num_days,
    num_infected_seeds, population_file, start_date, num_scenarios, num_runs, num_parallel_workers):

    event_log_level = "Transmissions"
    run_simplified = "false"
    track_index_case = "false"

    # Sample values for infectiousness and contact heterogeneity dispersion
    sampler = qmc.LatinHypercube(d=2)
    sample = sampler.random(n=num_scenarios)

    l_bounds = [0.2, 0.2]
    u_bounds = [0.6, 0.6]
    sample = qmc.scale(sample, l_bounds, u_bounds)

    #plt.scatter([v[0] for v in sample], [v[1] for v in sample])
    #plt.show()

    # Baseline
    run_parallel(scenario_name="baseline", contact_distribution="Constant", contact_distribution_overdispersion=0,
                    disease_config_file=disease_config_file, event_log_level=event_log_level,
                    holidays_file=holidays_file, num_days=num_days, num_infected_seeds=num_infected_seeds,
                    population_file=population_file, run_simplified=run_simplified, start_date=start_date,
                    track_index_case=track_index_case, transmission_probability_distribution="Constant",
                    transmission_probability=mean_transmission_probability,
                    transmission_probability_distribution_overdispersion=0, num_runs=num_runs,
                    num_parallel_workers=num_parallel_workers)

    scenario_names = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P"]
    for s_i in range(num_scenarios):
        scenario_name = scenario_names[s_i]
        infectiousness_overdispersion = sample[s_i][0]
        contacts_overdispersion = sample[s_i][1]

        mean_transmission_probability_corrected = get_mean_non_truncated_gamma(mean_transmission_probability, infectiousness_overdispersion)
        run_parallel(scenario_name=scenario_name, contact_distribution="Gamma", contact_distribution_overdispersion=contacts_overdispersion,
                        disease_config_file=disease_config_file, event_log_level=event_log_level,
                        holidays_file=holidays_file, num_days=num_days,
                        num_infected_seeds=num_infected_seeds, population_file=population_file,
                        run_simplified=run_simplified, start_date=start_date, track_index_case=track_index_case,
                        transmission_probability_distribution="Gamma", transmission_probability=mean_transmission_probability_corrected,
                        transmission_probability_distribution_overdispersion=infectiousness_overdispersion, num_runs=num_runs,
                        num_parallel_workers=num_parallel_workers)

if __name__=="__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--disease_config_file", type=str, default="disease_covid19_lognorm.xml")
    parser.add_argument("--holidays_file", type=str, default="holidays_belgium_2019_2021.csv")
    parser.add_argument("--transmission_probability", type=float, default=0.08)

    parser.add_argument("--num_days", type=int, default=200)
    parser.add_argument("--num_infected_seeds", type=int, default=1)
    parser.add_argument("--population_file", type=str, default="pop_belgium11M_c500_teachers_censushh.csv")
    parser.add_argument("--start_date", type=str, default="2020-02-17")

    parser.add_argument("--num_scenarios", type=int, default=6)
    parser.add_argument("--num_runs", type=int, default=4)
    parser.add_argument("--num_parallel_workers", type=int, default=4)

    args = parser.parse_args()
    run_lhs(args.disease_config_file, args.holidays_file, args.transmission_probability,
                args.num_days, args.num_infected_seeds, args.population_file,
                args.start_date, args.num_scenarios, args.num_runs, args.num_parallel_workers)
