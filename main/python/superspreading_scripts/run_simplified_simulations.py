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

import argparse

from run_util import get_mean_non_truncated_gamma, run_parallel

def run_simplified_simulations(disease_config_file, population_file, infected_seed_ids,
    start_date, num_runs, num_parallel_workers):

    event_log_level = "Transmissions"
    run_simplified = "true"
    track_index_case = "true"

    holidays_file = "holidays_none.csv"
    num_days = 40

    mean_transmission_probabilities = [0.025, 0.05, 0.075, 0.1]

    for person_id in infected_seed_ids:
        for tp in mean_transmission_probabilities:
            postfix = "_pid_" + str(person_id) + "_tp_" + str(tp)

            # Baseline
            run_parallel(scenario_name="simplified_baseline" + postfix,
                            contact_distribution="Constant", contact_distribution_overdispersion=0,
                            disease_config_file=disease_config_file, event_log_level=event_log_level,
                            holidays_file=holidays_file, num_days=num_days, num_infected_seeds=1,
                            population_file=population_file, run_simplified=run_simplified, start_date=start_date,
                            track_index_case=track_index_case, transmission_probability_distribution="Constant",
                            transmission_probability=tp, transmission_probability_distribution_overdispersion=0,
                            num_runs=num_runs, infected_seed_id=person_id, num_parallel_workers=num_parallel_workers)

            # Vary overdispersion infectiousness
            scenario_names = ["simplified_infectiousness_overdispersion_1000" + postfix,
                                "simplified_infectiousness_overdispersion_100" + postfix,
                                "simplified_infectiousness_overdispersion_60" + postfix,
                                "simplified_infectiousness_overdispersion_40" + postfix,
                                "simplified_infectiousness_overdispersion_20" + postfix]
            overdispersion_parameters = [10, 1, 0.6, 0.4, 0.2]

            for i in range(len(scenario_names)):
                scenario_name = scenario_names[i]
                k = overdispersion_parameters[i]
                tp_corrected = get_mean_non_truncated_gamma(tp, k, num_parallel_workers)
                run_parallel(scenario_name=scenario_name,
                                contact_distribution="Constant", contact_distribution_overdispersion=0,
                                disease_config_file=disease_config_file, event_log_level=event_log_level,
                                holidays_file=holidays_file, num_days=num_days, num_infected_seeds=1,
                                population_file=population_file, run_simplified=run_simplified, start_date=start_date,
                                track_index_case=track_index_case, transmission_probability_distribution="Gamma",
                                transmission_probability=tp_corrected, transmission_probability_distribution_overdispersion=k,
                                num_runs=num_runs, infected_seed_id=person_id, num_parallel_workers=num_parallel_workers)

            # Vary overdispersion contacts
            scenario_names = ["simplified_contacts_overdispersion_1000" + postfix,
                                "simplified_contacts_overdispersion_100" + postfix,
                                "simplified_contacts_overdispersion_60" + postfix,
                                "simplified_contacts_overdispersion_40" + postfix,
                                "simplified_contacts_overdispersion_20" + postfix]
            overdispersion_parameters = [10, 1, 0.6, 0.4, 0.2]

            for i in range(len(scenario_names)):
                scenario_name = scenario_names[i]
                k = overdispersion_parameters[i]
                run_parallel(scenario_name=scenario_name,
                                contact_distribution="Gamma", contact_distribution_overdispersion=k,
                                disease_config_file=disease_config_file, event_log_level=event_log_level,
                                holidays_file=holidays_file, num_days=num_days, num_infected_seeds=1,
                                population_file=population_file, run_simplified=run_simplified, start_date=start_date,
                                track_index_case=track_index_case, transmission_probability_distribution="Constant",
                                transmission_probability=tp, transmission_probability_distribution_overdispersion=0,
                                num_runs=num_runs, infected_seed_id=person_id, num_parallel_workers=num_parallel_workers)


if __name__=="__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--disease_config_file", type=str, default="disease_covid19_lognorm_nocntreduction.xml")
    parser.add_argument("--population_file", type=str, default="pop_belgium11M_c500_teachers_censushh.csv")
    parser.add_argument("--infected_seed_ids", type=int, nargs="+", default=[2,11,6]) # Adult (46), child (5), elderly (75)
    parser.add_argument("--start_date", type=str, default="2020-02-17")
    parser.add_argument("--num_runs", type=int, default=4)
    parser.add_argument("--num_parallel_workers", type=int, default=4)

    args = parser.parse_args()

    run_simplified_simulations(args.disease_config_file, args.population_file,
                                args.infected_seed_ids, args.start_date,
                                args.num_runs, args.num_parallel_workers)
