import argparse

from run_util import get_mean_non_truncated_gamma, run_parallel

def run_index_case_only_simulations(disease_config_file, num_infected_seeds, population_file, start_date, num_runs, num_parallel_workers):
    event_log_level = "Transmissions"
    run_simplified = "false"
    track_index_case = "true"

    holidays_file = "holidays_none.csv"
    num_days = 40

    mean_transmission_probabilities = [0.025, 0.05, 0.075, 0.1]

    for tp in mean_transmission_probabilities:
        # Baseline
        run_parallel(scenario_name="ico_baseline_tp_" + str(tp), contact_distribution="Constant", contact_distribution_overdispersion=0,
                        disease_config_file=disease_config_file, event_log_level=event_log_level,
                        holidays_file=holidays_file, num_days=num_days, num_infected_seeds=num_infected_seeds,
                        population_file=population_file, run_simplified=run_simplified, start_date=start_date,
                        track_index_case=track_index_case, transmission_probability_distribution="Constant",
                        transmission_probability=tp, transmission_probability_distribution_overdispersion=0,
                        num_runs=num_runs, num_parallel_workers=num_parallel_workers)

        # Vary dispersion parameter of infectiousness distribution
        scenario_names = ["ico_infectiousness_overdispersion_1000",
                            "ico_infectiousness_overdispersion_100",
                            "ico_infectiousness_overdispersion_60",
                            "ico_infectiousness_overdispersion_40",
                            "ico_infectiousness_overdispersion_20"]
        overdispersion_parameters = [10, 1, 0.6, 0.4, 0.2]

        for i in range(len(scenario_names)):
            scenario_name = scenario_names[i]
            k = overdispersion_parameters[i]
            mean_transmission_probability_corrected = get_mean_non_truncated_gamma(tp, k, num_parallel_workers)
            run_parallel(scenario_name=scenario_name + "_tp_" + str(tp), contact_distribution="Constant", contact_distribution_overdispersion=0,
                            disease_config_file=disease_config_file, event_log_level=event_log_level,
                            holidays_file=holidays_file, num_days=num_days, num_infected_seeds=num_infected_seeds,
                            population_file=population_file, run_simplified=run_simplified, start_date=start_date,
                            track_index_case=track_index_case, transmission_probability_distribution="Gamma",
                            transmission_probability=mean_transmission_probability_corrected,
                            transmission_probability_distribution_overdispersion=k, num_runs=num_runs,
                            num_parallel_workers=num_parallel_workers)

        # Vary dispersion parameter of contacts distribution
        scenario_names = ["ico_contacts_overdispersion_1000",
                            "ico_contacts_overdispersion_100",
                            "ico_contacts_overdispersion_60",
                            "ico_contacts_overdispersion_40",
                            "ico_contacts_overdispersion_20"]

        overdispersion_parameters = [10, 1, 0.6, 0.4, 0.2]

        for i in range(len(scenario_names)):
            scenario_name = scenario_names[i]
            k = overdispersion_parameters[i]

            run_parallel(scenario_name=scenario_name + "_tp_" + str(tp), contact_distribution="Gamma",
                            contact_distribution_overdispersion=k, disease_config_file=disease_config_file,
                            event_log_level=event_log_level, holidays_file=holidays_file, num_days=num_days,
                            num_infected_seeds=num_infected_seeds, population_file=population_file,
                            run_simplified=run_simplified, start_date=start_date, track_index_case=track_index_case,
                            transmission_probability_distribution="Constant", transmission_probability=tp,
                            transmission_probability_distribution_overdispersion=0, num_runs=num_runs,
                            num_parallel_workers=num_parallel_workers)


if __name__=="__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--disease_config_file", type=str, default="disease_covid19_lognorm.xml")
    parser.add_argument("--num_infected_seeds", type=int, default=1)
    parser.add_argument("--population_file", type=str, default="pop_belgium11M_c500_teachers_censushh.csv")
    parser.add_argument("--start_date", type=str, default="2020-02-17")
    parser.add_argument("--num_runs", type=int, default=4)
    parser.add_argument("--num_parallel_workers", type=int, default=4)

    args = parser.parse_args()

    run_index_case_only_simulations(args.disease_config_file, args.num_infected_seeds,
                                        args.population_file, args.start_date, args.num_runs,
                                        args.num_parallel_workers)
