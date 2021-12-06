import matplotlib.pyplot as plt
import numpy as np
import os

from collections import Counter
from scipy.stats import nbinom, probplot

def plot_ar(output_dir, fig_name, display_scenario_names, all_total_cases, num_days, pop_size, extinction_threshold=0, y_min=-0.05, y_max=1.05):
    all_total_cases_prop_pop = []
    means = []
    for scenario in all_total_cases:
        ars = [total_cases / pop_size for total_cases in scenario if total_cases >= extinction_threshold]
        all_total_cases_prop_pop.append(ars)
        means.append(np.mean(ars))

    plt.scatter(range(1, len(all_total_cases_prop_pop) + 1), means, marker="o")
    plt.boxplot(all_total_cases_prop_pop, labels=display_scenario_names)

    plt.ylabel("AR (after {} days)".format(num_days))
    plt.ylim(y_min, y_max)

    save_figure(output_dir, fig_name)

def plot_cumulative_cases_per_day(output_dir, fig_name, scenario_name, cases_by_day, num_days, y_max=None, events={}):
    for run in cases_by_day:
        total_cases = 0
        cumulative_cases_by_day = []
        for day in range(num_days):
            total_cases += run[day] if day in run else 0
            cumulative_cases_by_day.append(total_cases)
        plt.plot(range(num_days), cumulative_cases_by_day)

    for event_name, event_day in events.items():
        plt.axvline(event_day, color="lightgrey")
        plt.text(event_day + 1, y_max - (y_max / 5), event_name, rotation=90, color="lightgrey")

    plt.xlabel("Simulation day")
    plt.ylabel("Cumulative cases")

    if y_max is not None:
        plt.ylim(0, y_max)

    save_figure(output_dir, fig_name + "_" + scenario_name)

def plot_day_of_last_infection(output_dir, fig_name, display_scenario_names, all_days_last_infection, num_days):
    plt.boxplot(all_days_last_infection, labels=display_scenario_names)

    #plt.xticks(range(1, len(all_days_last_infection) + 1), display_scenario_names)

    plt.ylabel("Day with last infection")
    plt.ylim(-0.5, num_days + 1)

    save_figure(output_dir, fig_name, extension="png")

def plot_day_of_peak(output_dir, fig_name, display_scenario_names, all_cases_per_day, num_days, extinction_threshold=0):
    all_peak_days = []
    for scenario in all_cases_per_day:
        peak_days = []
        for run in scenario:
            if sum(run.values()) >= extinction_threshold:
                peak_days.append(max(run, key=lambda day: run[day]))
        all_peak_days.append(peak_days)
    plt.boxplot(all_peak_days, labels=display_scenario_names)

    plt.ylabel("Day of peak")
    plt.ylim(-0.5, num_days + 1)

    save_figure(output_dir, fig_name)

def plot_effective_r_by_day(output_dir, fig_name, scenario_name, rt_by_day, num_days, y_max=None, events={}):
    mean = []
    lower = []
    upper = []

    for day in range(num_days):
        r_on_day = [run[day] for run in rt_by_day]
        mean.append(np.nanmean(r_on_day))
        lower.append(np.percentile(r_on_day, 2.5))
        upper.append(np.percentile(r_on_day, 97.5))

    for event_name, event_day in events.items():
        plt.axvline(event_day, color="lightgrey")
        plt.text(event_day + 1, y_max - (y_max / 5), event_name, rotation=90, color="lightgrey")

    plt.plot(range(num_days), mean)
    plt.fill_between(range(num_days), lower, upper, color="lightgrey")

    plt.plot(range(num_days), [1] * num_days, color="orange") # Reference line at Rt = 1

    plt.xlabel("Simulation day")
    plt.xlim(-0.5, num_days)

    plt.ylabel("Rt")
    if y_max is not None:
        plt.ylim(0, y_max)

    save_figure(output_dir, fig_name + "_" + scenario_name, extension="png", dpi=100)

def plot_extinction_probabilities(output_dir, fig_name, display_scenario_names, all_total_cases, extinction_threshold, ylabel="Extinction probability"):
    extinction_probabilities = []
    for scenario in all_total_cases:
        extinction_probabilities.append(len([x for x in scenario if x < extinction_threshold]) / len(scenario))

    plt.bar(range(len(all_total_cases)), extinction_probabilities)
    plt.xticks(range(len(display_scenario_names)), display_scenario_names)
    plt.ylabel(ylabel + " (threshold = {} cases)".format(extinction_threshold))
    plt.ylim(0, 1.1)

    save_figure(output_dir, fig_name)

def plot_final_size_frequencies(output_dir, fig_name, display_scenario_names, all_total_cases, xlabel):
    """
        Plot frequency with which final sizes occur.
        Visualisation for extinction threshold.
    """

    plt.hist(all_total_cases, histtype="bar", stacked=True)

    plt.xlabel(xlabel)
    plt.ylabel("Frequency")
    plt.legend(display_scenario_names)

    save_figure(output_dir, fig_name)

def plot_herd_immunity_threshold(output_dir, fig_name, display_scenario_names, all_hits, show_day=False, num_days=200):
    hits = []
    for scenario in all_hits:
        hits.append([x for x in scenario if not np.isnan(x)])

    plt.boxplot(hits, labels=display_scenario_names)

    if show_day:
        plt.ylabel("Day on which Rt >= 1 for the last time")
        plt.ylim(0, num_days)
    else:
        plt.ylabel("Herd immunity threshold")
        plt.ylim(0, 1)

    save_figure(output_dir, fig_name)

def plot_new_cases_per_day(output_dir, fig_name, scenario_name, cases_by_day, num_days, y_max=None, events={}):
    for run in cases_by_day:
        plt.plot(range(num_days), [run[day] if day in run else 0 for day in range(num_days)])

    for event_name, event_day in events.items():
        plt.axvline(event_day, color="lightgrey")
        plt.text(event_day + 1, y_max - (y_max / 5), event_name, rotation=90, color="lightgrey")
    plt.xlabel("Simulation day")
    plt.ylabel("New cases")
    if y_max is not None:
        plt.ylim(0, y_max)

    save_figure(output_dir, fig_name + "_" + scenario_name)

def plot_num_cases_over_period(output_dir, fig_name, display_scenario_names, start_day, end_day, all_cases_over_period):
    #plt.boxplot(all_cases_over_period, labels=display_scenario_names)
    plt.violinplot(all_cases_over_period)
    plt.scatter(range(1, len(all_cases_over_period) + 1), [np.mean(scenario) for scenario in all_cases_over_period])

    plt.xticks(range(1, len(all_cases_over_period) + 1), display_scenario_names)
    plt.ylabel("Number of cases day {}-{}".format(start_day, end_day))

    save_figure(output_dir, fig_name, extension="png")

def plot_offspring_distributions(output_dir, fig_name, scenario_name, secondary_cases_by_tp):
    for tp in secondary_cases_by_tp:
        num_runs = len(secondary_cases_by_tp[tp])
        freq = Counter(secondary_cases_by_tp[tp])
        num_cases_sorted = list(freq.keys())
        num_cases_sorted.sort()

        plt.plot(num_cases_sorted, [freq[num] / num_runs for num in num_cases_sorted], marker="o")

    plt.xlabel("Number of secondary cases")
    plt.ylabel("Frequency")
    plt.ylim(-0.05, 1.05)
    plt.legend(["{:.3f}".format(tp) for tp in secondary_cases_by_tp], title="E(Transmission probability)")

    save_figure(output_dir, fig_name + "_" + scenario_name)

def plot_p80s(output_dir, fig_name, display_scenario_names, p80s):
    # Remove NaNs
    p80s = [[p80 for p80 in scenario_result if not np.isnan(p80)] for scenario_result in p80s]
    # Create boxplots
    plt.boxplot(p80s, labels=display_scenario_names)
    plt.ylabel("P80")
    plt.ylim(0, 1.1)
    save_figure(output_dir, fig_name, extension="png", dpi=100)

def plot_peak_sizes(output_dir, fig_name, display_scenario_names, all_cases_by_day, extinction_threshold=0, ymin=None, ymax=None):
    all_peak_sizes = []
    for scenario in all_cases_by_day:
        all_peak_sizes.append([max(run.values()) for run in scenario if sum(run.values()) >= extinction_threshold])
    plt.boxplot(all_peak_sizes, labels=display_scenario_names)
    if (ymin is not None) and (ymax is not None):
        plt.ylim(ymin, ymax)
    plt.ylabel("Peak size")

    save_figure(output_dir, fig_name)

def plot_qq(output_dir, fig_name, scenario_name, k, secondary_cases_by_tp):
    """
        Create QQ-plots to compare offspring distribution
        to negative binomial distribution
    """
    tps_sorted = list(secondary_cases_by_tp.keys())
    tps_sorted.sort()

    for tp in tps_sorted:
        res = probplot(secondary_cases_by_tp[tp], dist=nbinom, sparams=(k, k / (np.mean(secondary_cases_by_tp[tp]) + k)), fit=False, plot=plt)
        plt.xlabel("Negative binomial distribution quantiles")
        plt.ylabel("Simulations results qunatiles") # FIXME Find better axis labels
        plt.title("")

        save_figure(output_dir, fig_name + "_" + scenario_name + "_tp_" + "{:.3f}".format(tp))

def plot_secondary_cases_per_index_case(output_dir, fig_name, display_scenario_names, all_secondary_cases_by_tp, exclude_extinction=False):
    i = 0
    for scenario in all_secondary_cases_by_tp:

        if exclude_extinction:
            scenario_output = {}
            for tp, secondary_cases in scenario.items():
                secondary_cases_exclude_extinction = [num_cases for num_cases in scenario[tp] if num_cases > 0]
                if len(secondary_cases_exclude_extinction) > 0:
                    scenario_output[tp] = secondary_cases_exclude_extinction
        else:
            scenario_output = scenario

        tp_sorted = list(scenario_output.keys())
        tp_sorted.sort()


        lower = [np.percentile(scenario_output[tp], 2.5) for tp in tp_sorted]
        upper = [np.percentile(scenario_output[tp], 97.5) for tp in tp_sorted]

        plt.plot(tp_sorted, [np.mean(scenario_output[tp]) for tp in tp_sorted], color="C"+str(i), label=display_scenario_names[i])
        plt.fill_between(tp_sorted, lower, upper, facecolor="C"+str(i), alpha=0.3)

        i += 1

    plt.legend()
    plt.xlabel("Mean transmission probability")
    plt.xlim(0.02, 0.105)
    plt.ylabel("Number of secondary cases per index case")
    save_figure(output_dir, fig_name, extension="png")

def plot_secondary_cases_distribution(output_dir, fig_name, display_scenario_names, all_secondary_cases_frequencies):
    for scenario in all_secondary_cases_frequencies:
        all_frequencies = {}
        for run in scenario:
            for num_secondary_cases, freq in run.items():
                if num_secondary_cases in all_frequencies:
                    all_frequencies[num_secondary_cases] += freq
                else:
                    all_frequencies[num_secondary_cases] = freq

        num_cases_sorted = list(all_frequencies.keys())
        num_cases_sorted.sort()

        plt.plot(num_cases_sorted, [all_frequencies[num] for num in num_cases_sorted])
    plt.xlabel("Number of secondary cases")
    plt.xlim(-5, 105)

    plt.ylabel("Frequency")
    plt.yscale("log")

    plt.legend(display_scenario_names)

    save_figure(output_dir, fig_name)


def plot_transmissions_by_location(output_dir, fig_name, scenario_name, transmissions_by_location):
    transmissions_by_location_dict = {
        "Household": [],
        "K12School": [],
        "College": [],
        "Workplace": [],
        "PrimaryCommunity": [],
        "SecondaryCommunity": [],
    }

    for run in transmissions_by_location:
        total_transmissions = sum(run.values())

        for location, value in run.items():
            if location in transmissions_by_location_dict:
                if total_transmissions > 0:
                    transmissions_by_location_dict[location].append(value / total_transmissions)
                else:
                    transmissions_by_location_dict[location].append(0)
            else:
                if total_transmissions > 0:
                    transmissions_by_location_dict[location] = [value / total_transmissions]
                else:
                    transmissions_by_location_dict[location].append(0)

    locations = list(transmissions_by_location_dict.keys())
    plt.boxplot([transmissions_by_location_dict[loc] for loc in locations], labels=locations)

    plt.xticks(rotation=45)

    plt.ylabel("Fraction of infections")
    plt.ylim(-0.05, 1)

    save_figure(output_dir, fig_name + "_" + scenario_name)

def save_figure(output_dir, figure_name, extension="eps", dpi=200):
    if not os.path.exists(os.path.join(output_dir, "fig")):
        os.mkdir(os.path.join(output_dir, "fig"))

    plt.savefig(os.path.join(output_dir, "fig", figure_name + "." + extension), format=extension, bbox_inches='tight', dpi=dpi)
    plt.clf()
