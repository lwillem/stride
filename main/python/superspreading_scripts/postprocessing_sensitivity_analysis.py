import argparse
import matplotlib.pyplot as plt
import multiprocessing
import numpy as np
import os

from plots import save_figure
from postprocessing_util import get_day_of_last_infection, get_experiment_ids, get_herd_immunity_threshold, get_output_per_day, get_summary_output, get_total_cases


def plot_means(output_dir, fig_name, all_means, display_scenario_names, ylabel, legend, y_min=0, y_max=1):
    for scenario in all_means:
        plt.plot(range(len(scenario)), scenario)
    plt.xticks(range(len(display_scenario_names)), display_scenario_names)
    plt.ylim(y_min, y_max)
    plt.ylabel(ylabel)
    plt.legend(legend)
    save_figure(output_dir, fig_name)

def main(output_dir, sa_scenario_name, num_seeds_values, tp_values, legend, num_parallel_workers):
    if sa_scenario_name == "sa_num_seeds":
        legend = ["Number of infected seed = {}".format(value) for value in legend]
    elif sa_scenario_name == "sa_tp":
        legend = ["Mean transmission probability = {}".format(value) for value in legend]
    extinction_threshold = 50

    overdispersion_scenario_names = ["infectiousness_overdispersion", "contacts_overdispersion"]
    overdispersion_parameters = ["1000", "100", "60", "40", "20"]
    for overdispersion_scenario in overdispersion_scenario_names:
        alpha = r"$\alpha_{i}$"
        if overdispersion_scenario == "contacts_overdispersion":
            alpha = r"$\alpha_{c}$"
        display_scenario_names = ["Baseline", alpha + " = 10", alpha + " = 1", alpha + " = 0.6", alpha + " = 0.4", alpha + " = 0.2"]

        all_mean_ars = []
        all_mean_hits = []
        all_mean_hit_days = []

        all_cases_per_day = []

        all_extinction_probabilities = []

        for num_seeds in num_seeds_values:
            for tp in tp_values:
                output_prefix = "sa_num_seeds_{}_tp_{}".format(num_seeds, tp)
                output_dir_full = os.path.join(output_dir, output_prefix)
                scenario_names = [output_prefix + "baseline"] + [output_prefix + overdispersion_scenario + "_" + alpha for alpha in overdispersion_parameters]

                mean_ar = []
                mean_hit = []
                mean_hit_day = []

                num_cases_per_day = []

                extinction_probabilities = []

                for scenario_name in scenario_names:
                    experiment_ids = get_experiment_ids(output_dir_full, scenario_name)
                    with multiprocessing.Pool(processes=num_parallel_workers) as pool:
                        output_summary = pool.starmap(get_summary_output, [(output_dir_full, scenario_name, exp_id) for exp_id in experiment_ids])
                        output_per_day = pool.starmap(get_output_per_day, [(output_dir_full, scenario_name, exp_id) for exp_id in experiment_ids])

                        num_days = output_summary[0]["num_days"]
                        population_size = output_summary[0]["population_size"]

                        # Sort total cases from high to low & print
                        # Used to determine extinction threshold
                        total_cases = [get_total_cases(run["cases_per_day"], num_days) for run in output_per_day]
                        total_cases.sort(reverse=True)
                        print(total_cases)

                        extinction_probabilities.append(len([cases for cases in total_cases if cases >= extinction_threshold]) / len(total_cases))

                        mean_ar.append(np.mean([cases / population_size for cases in total_cases if cases >= extinction_threshold]))
                        mean_hit.append(np.mean([get_herd_immunity_threshold(run["rt_by_day"], run["cases_per_day"], num_days, population_size) for run in output_per_day if sum(run["cases_per_day"].values()) >= extinction_threshold]))
                        mean_hit_day.append(np.mean([get_herd_immunity_threshold(run["rt_by_day"], run["cases_per_day"], num_days, population_size, get_day=True) for run in output_per_day if sum(run["cases_per_day"].values()) >= extinction_threshold]))

                        num_cases_per_day.append([run["cases_per_day"] for run in output_per_day])
                all_mean_ars.append(mean_ar)
                all_mean_hits.append(mean_hit)
                all_mean_hit_days.append(mean_hit_day)
                all_cases_per_day.append(num_cases_per_day)
                all_extinction_probabilities.append(extinction_probabilities)

        plot_means(output_dir, "ar_" + sa_scenario_name + "_" + overdispersion_scenario, all_mean_ars, display_scenario_names, "mean AR after {} days".format(num_days), legend, y_min=0.5, y_max=1)
        plot_means(output_dir, "hit_" + sa_scenario_name + "_" + overdispersion_scenario, all_mean_hits, display_scenario_names, "mean herd immunity threshold", legend)
        plot_means(output_dir, "hit_day_" + sa_scenario_name + "_" + overdispersion_scenario, all_mean_hit_days, display_scenario_names, "mean day on which herd immunity threshold is reached", legend, y_min=0, y_max=200)
        plot_means(output_dir, "extinction_probabilities_" + sa_scenario_name + "_" + overdispersion_scenario, all_extinction_probabilities, display_scenario_names, "extinction probability", legend, y_min=0, y_max=1.05)

        all_mean_peak_sizes = []
        all_mean_peak_days = []
        all_mean_days_of_last_infection = []
        for sa_scenario in all_cases_per_day:
            mean_peak_size = []
            mean_peak_day = []
            mean_day_of_last_infection = []
            for scenario in sa_scenario:
                mean_peak_size.append(np.mean([max(run.values()) for run in scenario if sum(run.values()) >= extinction_threshold]))
                mean_peak_day.append(np.mean([max(run, key=lambda day: run[day]) for run in scenario if sum(run.values()) >= extinction_threshold]))
                mean_day_of_last_infection.append(np.mean([get_day_of_last_infection(run, num_days) for run in scenario if sum(run.values()) >= extinction_threshold]))
            all_mean_peak_sizes.append(mean_peak_size)
            all_mean_peak_days.append(mean_peak_day)
            all_mean_days_of_last_infection.append(mean_day_of_last_infection)

        plot_means(output_dir, "peak_size_" + sa_scenario_name + "_" + overdispersion_scenario, all_mean_peak_sizes, display_scenario_names, "mean peak size", legend, y_min=200000, y_max=800000)
        plot_means(output_dir, "peak_day_" + sa_scenario_name + "_" + overdispersion_scenario, all_mean_peak_days, display_scenario_names, "mean day of peak", legend, y_min=0, y_max=200)
        plot_means(output_dir, "day_of_last_infection_" + sa_scenario_name + "_" + overdispersion_scenario, all_mean_days_of_last_infection, display_scenario_names, "mean day of last infection", legend, y_min=0, y_max=200)

if __name__=="__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("output_dir", type=str)
    parser.add_argument("--sa_scenario_name", type=str, default="sa")
    parser.add_argument("--num_seeds_values", type=int, nargs="+", default=[1])
    parser.add_argument("--tp_values", type=float, nargs="+", default=[0.08])
    parser.add_argument("--legend", type=str, nargs="+", default=[])
    parser.add_argument("--num_parallel_workers", type=int, default=4)

    args = parser.parse_args()
    main(args.output_dir, args.sa_scenario_name, args.num_seeds_values, args.tp_values, args.legend, args.num_parallel_workers)
