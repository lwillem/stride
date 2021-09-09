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
Evolution of number of cases
Evolution of Rt
Theoretical description comparison
Index case only
R0
Distribution of number of secondary cases
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

def plot_new_cases_per_day(output_dir, fig_name, display_scenario_names, cases_by_day, num_days):
    for run in cases_by_day:
        plt.plot(range(num_days), [run[day] if day in run else 0 for day in range(num_days)])
    plt.show()

'''
def plot_evolution(output_dir, fig_name, scenario_name, new_cases_per_day, num_days, y_max):
    for run in new_cases_per_day:
        plt.plot(range(num_days), [run[day] for day in range(num_days)])

    plt.xlabel("Simulation day")
    plt.ylabel("New cases")
    plt.ylim(0, y_max)

    save_figure(output_dir, fig_name + "_" + scenario_name)
'''

def save_figure(output_dir, figure_name, extension="eps", dpi=1000):
    if not os.path.exists(os.path.join(output_dir, "fig")):
        os.mkdir(os.path.join(output_dir, "fig"))

    plt.savefig(os.path.join(output_dir, "fig", figure_name + "." + extension), format=extension, bbox_inches='tight', dpi=dpi)
    plt.clf()
