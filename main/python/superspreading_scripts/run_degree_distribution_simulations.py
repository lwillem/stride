#################################################################################
# Run simulations to get degree distribution.                                   #
#################################################################################

import argparse

from run_util import run_parallel

def run_degree_distribution_simulations(output_prefix, disease_config_file, num_days, population_file, start_date, num_parallel_workers):
    event_log_level = "All"
    run_simplified = "false"
    track_index_case = "false"

    holidays_file = "holidays_none.csv"
    mean_transmission_probability = 0.0 # No transmissions simulated for these experiments
    num_infected_seeds = 0

    num_runs = 1

    # Baseline
    run_parallel(scenario_name=output_prefix + "dd_baseline", contact_distribution="Constant", contact_distribution_overdispersion=0,
                    disease_config_file=disease_config_file, event_log_level=event_log_level,
                    holidays_file=holidays_file, num_days=num_days, num_infected_seeds=num_infected_seeds,
                    population_file=population_file, run_simplified=run_simplified, start_date=start_date,
                    track_index_case=track_index_case, transmission_probability_distribution="Constant",
                    transmission_probability=mean_transmission_probability,
                    transmission_probability_distribution_overdispersion=0, num_runs=num_runs,
                    num_parallel_workers=num_parallel_workers, summarize="DegreeDistribution")

    # Vary dispersion parameter of contacts distribution
    scenario_names = ["dd_contacts_overdispersion_1000",
                        "dd_contacts_overdispersion_100",
                        "dd_contacts_overdispersion_60",
                        "dd_contacts_overdispersion_40",
                        "dd_contacts_overdispersion_20"]

    overdispersion_parameters = [10, 1, 0.6, 0.4, 0.2]

    for i in range(len(scenario_names)):
        scenario_name = scenario_names[i]
        k = overdispersion_parameters[i]
        run_parallel(scenario_name=output_prefix + scenario_name, contact_distribution="Gamma",
                        contact_distribution_overdispersion=k, disease_config_file=disease_config_file,
                        event_log_level=event_log_level, holidays_file=holidays_file, num_days=num_days,
                        num_infected_seeds=num_infected_seeds, population_file=population_file,
                        run_simplified=run_simplified, start_date=start_date, track_index_case=track_index_case,
                        transmission_probability_distribution="Constant",
                        transmission_probability=mean_transmission_probability,
                        transmission_probability_distribution_overdispersion=0, num_runs=num_runs,
                        num_parallel_workers=num_parallel_workers, summarize="DegreeDistribution")


if __name__=="__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--output_prefix", type=str, default="")
    parser.add_argument("--disease_config_file", type=str, default="disease_covid19_lognorm.xml")
    parser.add_argument("--num_days", type=int, default=7)
    parser.add_argument("--population_file", type=str, default="pop_belgium11M_c500_teachers_censushh.csv")
    parser.add_argument("--start_date", type=str, default="2020-02-17")
    parser.add_argument("--num_parallel_workers", type=int, default=4)

    args = parser.parse_args()

    run_degree_distribution_simulations(args.output_prefix, args.disease_config_file, args.num_days, args.population_file, args.start_date,
                                            args.num_parallel_workers)


# TODO degree distribution when social distancing is practiced
