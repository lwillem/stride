'''
Degree distribution
Distribution of number of secondary cases
Extinction probability (heatmap)
Resurgence probability (heatmap)
Day of peak
Day of last infection
Herd immunity threshold
Total number of cases / AR
Over entire epidemic
Before lockdown /  during lockdown / after lockdown
Evolution of cumulative cases
Evolution of Rt
Theoretical description comparison
Index case only
R0
Distribution of number of secondary cases
'''

'''
def plot_ar(output_dir, fig_name, display_scenario_names, all_total_cases, num_days, pop_size, extinction_threshold=0, violin_plot=False):
    all_total_cases_prop_pop = []
    for scenario in all_total_cases:
        all_total_cases_prop_pop.append([total_cases / pop_size for total_cases in scenario if total_cases >= extinction_threshold])

    if violin_plot:
        plt.violinplot(all_total_cases_prop_pop)
        plt.xticks(range(1, len(all_total_cases_prop_pop) + 1), display_scenario_names, rotation=25)
    else:
        plt.boxplot(all_total_cases_prop_pop, labels=display_scenario_names)
        plt.xticks(rotation=25)

    plt.ylabel("AR (after {} days)".format(num_days))

    save_figure(output_dir, fig_name, extension="png")

def plot_day_last_infection(output_dir, fig_name, display_scenario_names, all_last_days_with_infections, violin_plot=False):
    if violin_plot:
        plt.violinplot(all_last_days_with_infections)
        plt.xticks(range(1, len(all_last_days_with_infections) + 1), display_scenario_names, rotation=25)
    else:
        plt.boxplot(all_last_days_with_infections, labels=display_scenario_names)
        plt.xticks(rotation=25)

    plt.plot(range(1, len(all_last_days_with_infections) + 1), [np.mean(x) for x in all_last_days_with_infections], linestyle="None", marker="o")

    plt.ylabel("Last day with new infections")
    save_figure(output_dir, fig_name, extension="png")
'''

import matplotlib.pyplot as plt
import numpy as np
import os

def plot_p80s(output_dir, fig_name, display_scenario_names, p80s):
    # Remove NaNs
    p80s = [[p80 for p80 in scenario_result if not np.isnan(p80)] for scenario_result in p80s]
    # Create boxplots
    plt.boxplot(p80s, labels=display_scenario_names)
    plt.ylabel("P80")
    plt.ylim(0, 0.5)
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

def plot_cumulative_cases_per_day(output_dir, fig_name, scenario_name, cases_by_day, num_days):
    for run in cases_by_day:
        total_cases = 0
        cumulative_cases_by_day = []
        for day in range(num_days):
            total_cases += run[day] if day in run else 0
            cumulative_cases_by_day.append(total_cases)
        plt.plot(range(num_days), cumulative_cases_by_day)

    plt.show()

def plot_new_cases_per_day(output_dir, fig_name, scenario_name, cases_by_day, num_days):
    for run in cases_by_day:
        plt.plot(range(num_days), [run[day] if day in run else 0 for day in range(num_days)])

    plt.xlabel("Simulation day")
    plt.ylabel("New cases")
    # TODO y_max ?

    save_figure(output_dir, fig_name + "_" + scenario_name, extension="png", dpi=100)

def save_figure(output_dir, figure_name, extension="eps", dpi=1000):
    if not os.path.exists(os.path.join(output_dir, "fig")):
        os.mkdir(os.path.join(output_dir, "fig"))

    plt.savefig(os.path.join(output_dir, "fig", figure_name + "." + extension), format=extension, bbox_inches='tight', dpi=dpi)
    plt.clf()
