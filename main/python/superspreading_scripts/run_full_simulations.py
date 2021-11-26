import argparse

from run_util import create_config, get_mean_non_truncated_gamma, run_parallel

def run_full_simulations(output_prefix, disease_config_file, holidays_file,
                            mean_transmission_probability, num_days, num_infected_seeds,
                            population_file, start_date, num_runs, num_parallel_workers):
    event_log_level = "Transmissions"
    run_simplified = "false"
    track_index_case = "false"

    # Baseline
    run_parallel(scenario_name=output_prefix + "baseline", contact_distribution="Constant", contact_distribution_overdispersion=0,
                    disease_config_file=disease_config_file, event_log_level=event_log_level,
                    holidays_file=holidays_file, num_days=num_days, num_infected_seeds=num_infected_seeds,
                    population_file=population_file, run_simplified=run_simplified, start_date=start_date,
                    track_index_case=track_index_case, transmission_probability_distribution="Constant",
                    transmission_probability=mean_transmission_probability,
                    transmission_probability_distribution_overdispersion=0, num_runs=num_runs,
                    num_parallel_workers=num_parallel_workers)

    # Vary dispersion parameter of infectiousness distribution
    scenario_names = ["infectiousness_overdispersion_1000",
                        "infectiousness_overdispersion_100",
                        "infectiousness_overdispersion_60",
                        "infectiousness_overdispersion_40",
                        "infectiousness_overdispersion_20"]

    overdispersion_parameters = [10, 1, 0.6, 0.4, 0.2]
    for i in range(len(scenario_names)):
        scenario_name = scenario_names[i]
        k = overdispersion_parameters[i]
        mean_transmission_probability_corrected = get_mean_non_truncated_gamma(mean_transmission_probability, k)
        run_parallel(scenario_name=output_prefix+scenario_name, contact_distribution="Constant", contact_distribution_overdispersion=0,
                        disease_config_file=disease_config_file,
                        event_log_level=event_log_level, holidays_file=holidays_file, num_days=num_days,
                        num_infected_seeds=num_infected_seeds, population_file=population_file,
                        run_simplified=run_simplified, start_date=start_date, track_index_case=track_index_case,
                        transmission_probability_distribution="Gamma", transmission_probability=mean_transmission_probability_corrected,
                        transmission_probability_distribution_overdispersion=k, num_runs=num_runs,
                        num_parallel_workers=num_parallel_workers)

    # Vary dispersion parameter of contacts distribution
    scenario_names = ["contacts_overdispersion_1000",
                        "contacts_overdispersion_100",
                        "contacts_overdispersion_60",
                        "contacts_overdispersion_40",
                        "contacts_overdispersion_20"]

    overdispersion_parameters = [10, 1, 0.6, 0.4, 0.2]
    for i in range(len(scenario_names)):
        scenario_name = scenario_names[i]
        k = overdispersion_parameters[i]
        run_parallel(scenario_name=output_prefix+scenario_name, contact_distribution="Gamma",
                            contact_distribution_overdispersion=k, disease_config_file=disease_config_file,
                            event_log_level=event_log_level, holidays_file=holidays_file, num_days=num_days,
                            num_infected_seeds=num_infected_seeds, population_file=population_file,
                            run_simplified=run_simplified, start_date=start_date, track_index_case=track_index_case,
                            transmission_probability_distribution="Constant",
                            transmission_probability=mean_transmission_probability,
                            transmission_probability_distribution_overdispersion=0, num_runs=num_runs,
                            num_parallel_workers=num_parallel_workers)

if __name__=="__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--output_prefix", type=str, default="")
    parser.add_argument("--disease_config_file", type=str, default="disease_covid19_lognorm.xml")
    parser.add_argument("--holidays_file", type=str, default="holidays_belgium_2019_2021.csv")
    parser.add_argument("--transmission_probability", type=float, default=0.08)
    parser.add_argument("--num_days", type=int, default=200)
    parser.add_argument("--num_infected_seeds", type=int, default=1)
    parser.add_argument("--population_file", type=str, default="pop_belgium11M_c500_teachers_censushh.csv")
    parser.add_argument("--start_date", type=str, default="2020-02-17")
    parser.add_argument("--num_runs", type=int, default=4)
    parser.add_argument("--num_parallel_workers", type=int, default=4)

    args = parser.parse_args()

    run_full_simulations(output_prefix=args.output_prefix, disease_config_file=args.disease_config_file,
                            holidays_file=args.holidays_file, mean_transmission_probability=args.transmission_probability,
                            num_days=args.num_days, num_infected_seeds=args.num_infected_seeds,
                            population_file=args.population_file,
                            start_date=args.start_date, num_runs=args.num_runs,
                            num_parallel_workers=args.num_parallel_workers)

# TODO combinations of infectiousness + contacts overdispersion?
