import matplotlib.pyplot as plt
import numpy as np
import os

from collections import Counter
from scipy.stats import nbinom, probplot
from statsmodels.nonparametric.smoothers_lowess import lowess
from statsmodels.stats.proportion import proportion_confint

from postprocessing_util import get_cases_over_period

def plot_ar(output_dir, fig_name, display_scenario_names, all_total_cases, num_days, pop_size, color, extinction_threshold=0, y_min=-0.05, y_max=1.05):
    all_total_cases_prop_pop = []
    means = []
    for scenario in all_total_cases:
        ars = [total_cases / pop_size for total_cases in scenario if total_cases >= extinction_threshold]
        all_total_cases_prop_pop.append(ars)
        means.append(np.mean(ars))

    violin_parts = plt.violinplot(all_total_cases_prop_pop)
    for part_name, part_properties in violin_parts.items():
        if part_name == "bodies":
            for pc in part_properties:
                pc.set_color(color)
        else:
            part_properties.set_color(color)

    plt.scatter(range(1, len(all_total_cases_prop_pop) + 1), means, color="orange")

    plt.xticks(range(1, len(all_total_cases_prop_pop) + 1), display_scenario_names)
    plt.ylabel("AR (after {} days)".format(num_days))
    plt.ylim(y_min, y_max)

    save_figure(output_dir, fig_name, extension="png")

def plot_cumulative_cases_per_day(output_dir, fig_name, scenario_name, cases_by_day, num_days, color, y_max=None, events={}):
    for run in cases_by_day:
        total_cases = 0
        cumulative_cases_by_day = []
        for day in range(num_days):
            total_cases += run[day] if day in run else 0
            cumulative_cases_by_day.append(total_cases)
        plt.plot(range(num_days), cumulative_cases_by_day, color=color)

    for event_name, event_day in events.items():
        plt.axvline(event_day, color="lightgrey")
        plt.text(event_day + 1, y_max - (y_max / 5), event_name, rotation=90, color="lightgrey")

    plt.xlabel("Simulation day")
    plt.ylabel("Cumulative cases")

    if y_max is not None:
        plt.ylim(0, y_max)

    save_figure(output_dir, fig_name + "_" + scenario_name)


def plot_day_of_last_infection(output_dir, fig_name, display_scenario_names, all_days_last_infection, y_min, num_days, color):
    violin_parts = plt.violinplot(all_days_last_infection)
    for part_name, part_properties in violin_parts.items():
        if part_name == "bodies":
            for pc in part_properties:
                pc.set_color(color)
        else:
            part_properties.set_color(color)

    plt.scatter(range(1, len(all_days_last_infection) + 1), [np.mean(scenario) for scenario in all_days_last_infection], color="orange")

    plt.xticks(range(1, len(all_days_last_infection) + 1), display_scenario_names)

    plt.ylabel("Day with last infection")
    plt.ylim(y_min, num_days + 1)

    save_figure(output_dir, fig_name, extension="png")

def plot_day_of_peak(output_dir, fig_name, display_scenario_names, all_cases_per_day, num_days, color, extinction_threshold=0):
    all_peak_days = []
    for scenario in all_cases_per_day:
        peak_days = []
        for run in scenario:
            if sum(run.values()) >= extinction_threshold:
                peak_days.append(max(run, key=lambda day: run[day]))
        all_peak_days.append(peak_days)

    violin_parts = plt.violinplot(all_peak_days)
    for part_name, part_properties in violin_parts.items():
        if part_name == "bodies":
            for pc in part_properties:
                pc.set_color(color)
        else:
            part_properties.set_color(color)

    plt.scatter(range(1, len(all_peak_days) + 1), [np.mean(scenario) for scenario in all_peak_days], color="orange")

    plt.xticks(range(1, len(all_peak_days) + 1), display_scenario_names)
    plt.ylabel("Day of peak")
    plt.ylim(-0.5, num_days + 1)

    save_figure(output_dir, fig_name, extension="png")

def plot_effective_r_by_day(output_dir, fig_name, scenario_name, rt_by_day, num_days, color, y_max=None, events={}, smoothed=False):

    if smoothed:
        for run_i in range(len(rt_by_day)):
            rt_by_day[run_i] = lowess([rt_by_day[run_i][day] for day in range(num_days)], range(num_days), is_sorted=True, return_sorted=False)

    mean = []
    lower = []
    upper = []

    for day in range(num_days):
        r_on_day = [run[day] for run in rt_by_day]
        r_on_day = [r for r in r_on_day if not np.isnan(r)]
        if len(r_on_day) > 0:
            mean.append(np.mean(r_on_day))
            lower.append(np.percentile(r_on_day, 2.5))
            upper.append(np.percentile(r_on_day, 97.5))
        else:
            mean.append(np.nan)
            lower.append(np.nan)
            upper.append(np.nan)

    for event_name, event_day in events.items():
        plt.axvline(event_day, color="darkgrey")
        plt.text(event_day + 2, y_max - (y_max / 5), event_name, rotation=90, color="darkgrey")

    plt.plot(range(num_days), mean, color=color)
    plt.fill_between(range(num_days), lower, upper, color="lightgrey")

    plt.plot(range(num_days), [1] * num_days, color="orange") # Reference line at Rt = 1

    plt.xlabel("Simulation day")
    plt.xlim(-0.5, num_days)

    plt.ylabel("Rt")
    if y_max is not None:
        plt.ylim(0, y_max)

    save_figure(output_dir, fig_name + "_" + scenario_name, extension="png", dpi=100)

def plot_extinction_probabilities(output_dir, fig_name, display_scenario_names, all_total_cases, extinction_threshold, color, ylabel="Extinction probability"):
    extinction_probabilities = []
    confidence_intervals = []
    for scenario in all_total_cases:
        num_runs = len(scenario)
        num_extinctions = len([x for x in scenario if x < extinction_threshold])
        extinction_probabilities.append(num_extinctions / num_runs)
        conf = proportion_confint(count=num_extinctions, nobs=num_runs, alpha=0.05, method="beta") # Clopper-Pearson interval
        confidence_intervals.append(conf)

    #print(extinction_probabilities)

    i = 0
    horizontal_line_width = 0.3
    for interval in confidence_intervals:
        bottom = interval[0]
        top = interval[1]
        left = i - horizontal_line_width / 2
        right = i + horizontal_line_width / 2
        plt.plot([i,i], [top, bottom], color=color)
        plt.plot([left, right], [top, top], color=color)
        plt.plot([left, right], [bottom, bottom], color=color)
        i += 1

    plt.plot(range(len(all_total_cases)), extinction_probabilities, marker="o", linestyle="None", color=color)

    plt.xticks(range(len(display_scenario_names)), display_scenario_names)
    plt.ylabel(ylabel + " (threshold = {} cases)".format(extinction_threshold))
    plt.ylim(0, 1.1)

    save_figure(output_dir, fig_name)

def plot_resurgence_probabilities(output_dir, fig_name, display_scenario_names, all_cases_per_day, start_lockdown, end_lockdown, num_days, resurgence_threshold, color):
    resurgence_probabilities = []
    confidence_intervals = []

    for scenario in all_cases_per_day:
        runs_with_cases_during_lockdown = 0
        runs_above_resurgence_threshold = 0

        for run in scenario:
            num_cases_during_lockdown = get_cases_over_period(run, start_lockdown, end_lockdown)
            num_cases_after_lockdown = get_cases_over_period(run, end_lockdown, num_days)
            if not num_cases_during_lockdown == 0:
                runs_with_cases_during_lockdown += 1
                if num_cases_after_lockdown >= resurgence_threshold:
                    runs_above_resurgence_threshold += 1
        resurgence_probabilities.append(runs_above_resurgence_threshold / runs_with_cases_during_lockdown)
        conf = proportion_confint(count=runs_above_resurgence_threshold, nobs=runs_with_cases_during_lockdown, alpha=0.05, method="beta") # Clopper-Pearson interval
        confidence_intervals.append(conf)

    # Plot confidence intervals
    i = 0
    horizontal_line_width = 0.3
    for interval in confidence_intervals:
        bottom = interval[0]
        top = interval[1]
        left = i - horizontal_line_width / 2
        right = i + horizontal_line_width / 2

        plt.plot([i,i], [top, bottom], color=color)
        plt.plot([left, right], [top, top], color=color)
        plt.plot([left, right], [bottom, bottom], color=color)
        i += 1

    # Plot resurgence probabilities
    plt.plot(range(len(resurgence_probabilities)), resurgence_probabilities, marker="o", linestyle="None", color=color)
    plt.xticks(range(len(display_scenario_names)), display_scenario_names)

    plt.ylabel("Resurgence probability")
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

def plot_herd_immunity_threshold(output_dir, fig_name, display_scenario_names, all_hits, color, show_day=False, num_days=200, y_min=0, y_max=1):
    hits = []
    means = []
    for scenario in all_hits:
        hits_no_nan = [x for x in scenario if not np.isnan(x)]
        means.append(np.mean(hits_no_nan))
        hits.append(hits_no_nan)

    violin_parts = plt.violinplot(hits)
    for part_name, part_properties in violin_parts.items():
        if part_name == "bodies":
            for pc in part_properties:
                pc.set_color(color)
        else:
            part_properties.set_color(color)
    plt.scatter(range(1, len(hits) + 1), means, color="orange")
    plt.xticks(range(1, len(hits) + 1), display_scenario_names)

    if show_day:
        plt.ylabel("Day on which Rt >= 1 for the last time")
        plt.ylim(0, num_days)
    else:
        plt.ylabel("Herd immunity threshold")
        plt.ylim(y_min, y_max)

    save_figure(output_dir, fig_name, extension="png")

def plot_new_cases_per_day(output_dir, fig_name, scenario_name, cases_by_day, num_days, color, y_max=None, events={}):
    for run in cases_by_day:
        plt.plot(range(num_days), [run[day] if day in run else 0 for day in range(num_days)], color=color)

    for event_name, event_day in events.items():
        plt.axvline(event_day, color="lightgrey")
        plt.text(event_day + 1, y_max - (y_max / 5), event_name, rotation=90, color="lightgrey")
    plt.xlabel("Simulation day")
    plt.ylabel("New cases")
    if y_max is not None:
        plt.ylim(0, y_max)

    save_figure(output_dir, fig_name + "_" + scenario_name)

def plot_num_cases_over_period(output_dir, fig_name, display_scenario_names, start_day, end_day, all_cases_over_period, color, y_min=None, y_max=None, extinction_threshold=0):
    violin_parts = plt.violinplot([[x for x in scenario if x >= extinction_threshold] for scenario in all_cases_over_period])
    for part_name, part_properties in violin_parts.items():
        if part_name == "bodies":
            for pc in part_properties:
                pc.set_color(color)
        else:
            part_properties.set_color(color)

    plt.scatter(range(1, len(all_cases_over_period) + 1), [np.mean([x for x in scenario if x >= extinction_threshold]) for scenario in all_cases_over_period], color="orange")

    plt.xticks(range(1, len(all_cases_over_period) + 1), display_scenario_names)
    plt.ylabel("Number of cases day {}-{}".format(start_day, end_day))
    if y_min is not None and y_max is not None:
        plt.ylim(y_min, y_max)

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

def plot_p80s(output_dir, fig_name, display_scenario_names, p80s, color):
    # Remove NaNs
    p80s = [[p80 for p80 in scenario_result if not np.isnan(p80)] for scenario_result in p80s]
    means = [np.mean(scenario_result) for scenario_result in p80s]
    print(means)

    violin_parts = plt.violinplot(p80s)
    for part_name, part_properties in violin_parts.items():
        if part_name == "bodies":
            for pc in part_properties:
                pc.set_color(color)
        else:
            part_properties.set_color(color)

    plt.xticks(range(1, len(display_scenario_names) + 1), display_scenario_names)

    plt.scatter(range(1, len(p80s) + 1), means, color="orange")

    plt.ylabel("P80")
    plt.ylim(0, 1.1)
    save_figure(output_dir, fig_name, extension="png", dpi=100)

def plot_peak_sizes(output_dir, fig_name, display_scenario_names, all_cases_by_day, start_day, end_day, color, extinction_threshold=0, ymin=None, ymax=None):
    all_peak_sizes = []
    means = []
    for scenario in all_cases_by_day:
        peak_sizes = [max(list(run.values())[start_day:end_day]) for run in scenario if sum(list(run.values())[start_day:end_day]) >= extinction_threshold]
        all_peak_sizes.append(peak_sizes)
        means.append(np.mean(peak_sizes))

    violin_parts = plt.violinplot(all_peak_sizes)
    for part_name, part_properties in violin_parts.items():
        if part_name == "bodies":
            for pc in part_properties:
                pc.set_color(color)
        else:
            part_properties.set_color(color)
    plt.scatter(range(1, len(all_peak_sizes) + 1), means, color="orange")

    plt.xticks(range(1, len(all_peak_sizes) + 1), display_scenario_names)

    if (ymin is not None) and (ymax is not None):
        plt.ylim(ymin, ymax)
    plt.ylabel("Peak size")

    save_figure(output_dir, fig_name, extension="png")

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
        plt.ylabel("Simulations results quantiles") # FIXME Find better axis labels
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

def plot_transmissions_by_location(output_dir, fig_name, scenario_name, transmissions_by_location, color):
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
    bplot = plt.boxplot([transmissions_by_location_dict[loc] for loc in locations], labels=locations, patch_artist=True)
    for part_name, part_list in bplot.items():
        if part_name == "boxes":
            for box in part_list:
                box.set_edgecolor(color)
                box.set_facecolor(color)
                box.set_alpha(0.3)
        else:
            for part in part_list:
                part.set_color(color)
                part.set_markeredgecolor(color)

    plt.xticks(rotation=45)

    plt.ylabel("Fraction of infections")
    plt.ylim(-0.05, 1.05)

    save_figure(output_dir, fig_name + "_" + scenario_name, extension="png")

def plot_secondary_cases_per_index_case_histogram(output_dir, fig_name, secondary_cases_per_index_case, display_scenario_names):
    plt.hist(secondary_cases_per_index_case, histtype="bar", stacked=True)
    plt.xlabel("Number of secondary cases")
    plt.ylabel("Frequency")
    plt.legend(display_scenario_names)

    save_figure(output_dir, fig_name)

def plot_secondary_cases_per_index_case_means(output_dir, fig_name, secondary_cases_per_index_case, display_scenario_names, color):
    horizontal_line_width = 0.5

    means = [np.mean(scenario) for scenario in secondary_cases_per_index_case]
    print(means)

    violin_parts = plt.violinplot(secondary_cases_per_index_case)
    for part_name, part_properties in violin_parts.items():
        if part_name == "bodies":
            for pc in part_properties:
                pc.set_color(color)
        else:
            part_properties.set_color(color)
    plt.scatter(range(1, len(secondary_cases_per_index_case) + 1), means, color="orange")

    plt.xticks(range(1, len(secondary_cases_per_index_case) + 1), display_scenario_names)
    plt.ylabel("Secondary cases per index case")
    plt.ylim(0, 40)

    save_figure(output_dir, fig_name, extension="png")

def plot_comparison_means_theoretical_sims(output_dir, fig_name, scenario_name, transmission_probabilities, secondary_cases_by_tp, theoretical_means):
    """
        Plot theoretical estimate of mean number of secondary cases per index case
        VS mean and 95% interval of number of secondary cases per index case from simulations.
    """
    means = [np.mean(secondary_cases) for secondary_cases in secondary_cases_by_tp]
    lower = [np.percentile(secondary_cases, 2.5) for secondary_cases in secondary_cases_by_tp]
    upper = [np.percentile(secondary_cases, 97.5) for secondary_cases in secondary_cases_by_tp]

    plt.plot(transmission_probabilities, means, marker="o", label="Simulations")
    plt.fill_between(transmission_probabilities, lower, upper, color="lightgrey")

    plt.plot(transmission_probabilities, theoretical_means, linestyle="None", marker="^", label="Theoretical")

    plt.xlabel("Mean transmission probality")
    plt.ylabel("Secondary cases per index case")
    plt.ylim(-2, 100) # TODO parameter?
    plt.legend()

    save_figure(output_dir, fig_name + "_" + scenario_name, extension="png")

def plot_comparison_variance_theoretical_sims(output_dir, fig_name, scenario_name, transmission_probabilities, secondary_cases_by_tp, theoretical_variances):
    """
        Plot theoretical estimate of variance for number of secondary cases per index case
        VS variance for number of secondary cases per index case from simulations.
    """
    variances = [np.var(secondary_cases) for secondary_cases in secondary_cases_by_tp]

    plt.plot(transmission_probabilities, variances, linestyle="None", marker="o")
    plt.plot(transmission_probabilities, theoretical_variances, linestyle="None", marker="^")

    plt.legend(["Simulations", "Theoretical"])

    plt.xlabel("Mean transmission probality")
    plt.ylabel("Variance of secondary cases caused by index case")
    save_figure(output_dir, fig_name + "_" + scenario_name, extension="png")

def save_figure(output_dir, figure_name, extension="eps", dpi=200):
    if not os.path.exists(os.path.join(output_dir, "fig")):
        os.mkdir(os.path.join(output_dir, "fig"))

    plt.savefig(os.path.join(output_dir, "fig", figure_name + "." + extension), format=extension, bbox_inches='tight', dpi=dpi)
    plt.clf()
