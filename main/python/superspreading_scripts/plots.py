import matplotlib.pyplot as plt
import numpy as np
import os

def plot_ar(output_dir, fig_name, display_scenario_names, all_total_cases, num_days, pop_size, extinction_threshold=0):
    all_total_cases_prop_pop = []
    for scenario in all_total_cases:
        all_total_cases_prop_pop.append([total_cases / pop_size for total_cases in scenario if total_cases >= extinction_threshold])

    plt.violinplot(all_total_cases_prop_pop)
    plt.xticks(range(1, len(all_total_cases_prop_pop) + 1), display_scenario_names)

    plt.ylabel("AR (after {} days)".format(num_days))
    plt.ylim(0, 1)

    save_figure(output_dir, fig_name, extension="png", dpi=100)

def plot_cumulative_cases_per_day(output_dir, fig_name, scenario_name, cases_by_day, num_days):
    for run in cases_by_day:
        total_cases = 0
        cumulative_cases_by_day = []
        for day in range(num_days):
            total_cases += run[day] if day in run else 0
            cumulative_cases_by_day.append(total_cases)
        plt.plot(range(num_days), cumulative_cases_by_day)

    plt.xlabel("Simulation day")
    plt.ylabel("Cumulative cases")

    # TODO y_max?

    save_figure(output_dir, fig_name + "_" + scenario_name, extension="png", dpi=100)

def plot_herd_immunity_threshold(output_dir, fig_name, display_scenario_names, all_hits, show_day=False):
    hits = []
    for scenario in all_hits:
        hits.append([x for x in scenario if not np.isnan(x)])

    plt.violinplot(hits)
    plt.xticks(range(1, len(display_scenario_names + 1), display_scenario_names))

    if show_day:
        plt.ylabel("Day on which Rt >= 1 for the last time")
    else:
        plt.ylabel("Herd immunity threshold")

    save_figure(output_dir, fig_name, extension="png", dpi=100)

def plot_extinction_probabilities(output_dir, fig_name, display_scenario_names, all_total_cases, extinction_threshold):
    extinction_probabilities = []
    for scenario in all_total_cases:
        extinction_probabilities.append(len([x for x in scenario if x < extinction_threshold]) / len(scenario))

    plt.bar(range(len(all_total_cases)), extinction_probabilities)
    plt.xticks(range(len(display_scenario_names)), display_scenario_names)
    plt.ylabel("Extinction probability (threshold = {} cases)".format(extinction_threshold))
    plt.ylim(0, 1.1)

    save_figure(output_dir, fig_name, extension="png", dpi=100)

def plot_final_size_frequencies(output_dir, fig_name, display_scenario_names, all_total_cases, num_days):
    """
        Plot frequency with which final sizes occur.
        Visualisation for extinction threshold.
    """

    plt.hist(all_total_cases, histtype="bar", stacked=True)

    plt.xlabel("Outbreak size after {} days".format(num_days))
    plt.ylabel("Frequency")
    plt.legend(display_scenario_names)

    save_figure(output_dir, fig_name, extension="png", dpi=100)

def plot_day_last_infection(output_dir, fig_name, display_scenario_names, all_days_last_infection, num_days):
    plt.violinplot(all_days_last_infection)

    plt.xticks(range(1, len(all_days_last_infection) + 1), display_scenario_names)

    plt.ylabel("Day with last infection")
    plt.ylim(0, num_days + 1)

    save_figure(output_dir, fig_name, extension="png", dpi=100)

def plot_new_cases_per_day(output_dir, fig_name, scenario_name, cases_by_day, num_days):
    for run in cases_by_day:
        plt.plot(range(num_days), [run[day] if day in run else 0 for day in range(num_days)])

    plt.xlabel("Simulation day")
    plt.ylabel("New cases")
    # TODO y_max ?

    save_figure(output_dir, fig_name + "_" + scenario_name, extension="png", dpi=100)

def plot_p80s(output_dir, fig_name, display_scenario_names, p80s):
    # Remove NaNs
    p80s = [[p80 for p80 in scenario_result if not np.isnan(p80)] for scenario_result in p80s]
    # Create boxplots
    plt.boxplot(p80s, labels=display_scenario_names)
    plt.ylabel("P80")
    plt.ylim(0, 0.5)
    save_figure(output_dir, fig_name, extension="png", dpi=100)

def plot_effective_r_by_day(output_dir, fig_name, scenario_name, rt_by_day, num_days):
    mean = []
    lower = []
    upper = []

    for day in range(num_days):
        r_on_day = [run[day] for run in rt_by_day]
        mean.append(np.nanmean(r_on_day))
        lower.append(np.percentile(r_on_day, 2.5))
        upper.append(np.percentile(r_on_day, 97.5))

    plt.plot(range(num_days), mean)
    plt.fill_between(range(num_days), lower, upper, color="lightgrey")

    plt.plot(range(num_days), [1] * num_days, color="orange") # Reference line at Rt = 1

    plt.xlabel("Simulation day")
    plt.xlim(-0.5, num_days)

    plt.ylabel("Rt")
    # TODO ylim?

    save_figure(output_dir, fig_name + "_" + scenario_name, extension="png", dpi=100)

def save_figure(output_dir, figure_name, extension="eps", dpi=1000):
    if not os.path.exists(os.path.join(output_dir, "fig")):
        os.mkdir(os.path.join(output_dir, "fig"))

    plt.savefig(os.path.join(output_dir, "fig", figure_name + "." + extension), format=extension, bbox_inches='tight', dpi=dpi)
    plt.clf()
