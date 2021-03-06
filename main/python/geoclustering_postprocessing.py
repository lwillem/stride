import argparse
import csv
import datetime
import matplotlib.pyplot as plt
import multiprocessing
import os

def get_experiment_ids(output_dir, scenario_name):
    experiment_ids = []
    summary_file = os.path.join(output_dir, scenario_name, scenario_name + "_summary.csv")
    with open(summary_file) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            experiment_ids.append(int(row["exp_id"]))
    return experiment_ids


def get_new_cases_per_day(output_dir, scenario_name, exp_id, num_days):
    transmissions_file = os.path.join(output_dir, scenario_name, "exp" + "{:04}".format(exp_id), "contact_log.txt")

    days = {}
    for d in range(num_days):
        days[d] = 0
    with open(transmissions_file) as f:
        for line in f:
            line = line.split()
            tag = line[0]
            if tag == "[TRAN]":
                sim_day = int(line[6])
                days[sim_day] += 1

    return days

def get_attack_rate(output_dir, scenario_name, exp_id, num_days):
    pop_size = 0
    summary_file = os.path.join(output_dir, scenario_name, "exp" + "{:04}".format(exp_id), "summary.csv")
    with open(summary_file) as csvfile:
        reader = csv.DictReader(csvfile)
        summary = next(reader)
        pop_size = int(summary["population_size"])

    total_cases = 0
    transmissions_file = os.path.join(output_dir, scenario_name, "exp" + "{:04}".format(exp_id), "contact_log.txt")
    with open(transmissions_file) as f:
        for line in f:
            line = line.split()
            tag = line[0]
            if tag == "[TRAN]":
                sim_day = int(line[6])
                if sim_day < num_days:
                    total_cases += 1
    if pop_size > 0:
        return total_cases / pop_size

def get_num_non_compliers(output_dir, scenario_name, exp_id):
    num_non_compliers = 0
    log_file = os.path.join(output_dir, scenario_name, "exp" + "{:04}".format(exp_id), "contact_log.txt")
    with open(log_file) as f:
        for line in f:
            line = line.split()
            tag = line[0]
            if tag == "[NCOM]":
                num_non_compliers += 1
    return num_non_compliers


def get_num_non_compliers_by_age(output_dir, scenario_name, exp_id):
    max_age = 111
    by_age = {}
    for a in range(max_age):
        by_age[a] = 0
    log_file = os.path.join(output_dir, scenario_name, "exp" + "{:04}".format(exp_id), "contact_log.txt")
    with open(log_file) as f:
        for line in f:
            line = line.split()
            tag = line[0]
            if tag == "[NCOM]":
                nc_age = int(line[2])
                by_age[nc_age] += 1

    return by_age

def main(output_dir, scenario_names, scenario_display_names):
    if scenario_display_names == None:
        scenario_display_names = scenario_names
    # FIXME get from config file?
    num_days = 120
    base = datetime.datetime.strptime("2020-08-31", "%Y-%m-%d")
    dates = [base + datetime.timedelta(days=x) for x in range(num_days)]

    all_new_cases_per_day = {}
    all_cumulative_cases_per_day = {}
    all_attack_rates = []
    all_nums_non_compliers = []
    all_nc_by_age = []

    for name in scenario_names:
        experiment_ids = get_experiment_ids(output_dir, name)
        with multiprocessing.Pool(processes=8) as pool:
            new_cases_per_day = pool.starmap(get_new_cases_per_day, [(output_dir, name, exp_id, num_days) for exp_id in experiment_ids])
            cumulative_cases_per_day = []
            for sim in new_cases_per_day:
                cumulative_cases_sim = {}
                total_cases = 0
                for day in sim:
                    total_cases += sim[day]
                    cumulative_cases_sim[day] = total_cases
                cumulative_cases_per_day.append(cumulative_cases_sim)

            all_new_cases_per_day[name] = new_cases_per_day
            all_cumulative_cases_per_day[name] = cumulative_cases_per_day

            '''by_age = pool.starmap(get_num_non_compliers_by_age, [(output_dir, name, exp_id) for exp_id in experiment_ids])
            mean_by_age = []
            for age in range(111):
                nc_for_age = [r[age] for r in by_age]
                mean = sum(nc_for_age) / len(nc_for_age)
                #mean = [sum(r[age]) / len(r[age]) for r in by_age]
                mean_by_age.append(mean)
            all_nc_by_age.append(mean_by_age)'''

    """for name in scenario_names:

            attack_rates = pool.starmap(get_attack_rate, [(output_dir, name, exp_id, num_days) for exp_id in experiment_ids])
            all_attack_rates.append(attack_rates)
            num_non_compliers = pool.starmap(get_num_non_compliers, [(output_dir, name, exp_id) for exp_id in experiment_ids])
            all_nums_non_compliers.append(num_non_compliers)

    """

    '''for i in range(len(scenario_names)):
        plt.plot(all_nc_by_age[i])
    plt.xlabel("Age (in years)")
    plt.ylabel("Mean number of non-compliers")
    plt.legend(scenario_names)

    plt.show()'''

    """
    plt.savefig("non_compliers_by_age_comparison.png")
    plt.clf()
    """

    i = 0
    colors = ["green", "red", "blue", "yellow", "magenta"]
    plots = []
    for name in scenario_names:
        results = all_new_cases_per_day[name]
        for r in results:
            new_cases = []
            for d in range(num_days):
                new_cases.append(r[d])
            p, = plt.plot(dates, new_cases, color=colors[i])
        plots.append(p)
        i += 1

    plt.xlabel("Day")
    plt.xticks(dates[::5], [x.strftime("%d/%m") for x in dates[::5]], rotation=45)
    plt.ylabel("New infections")

    plt.legend(plots, scenario_display_names)
    plt.tight_layout()
    plt.show()


    i = 0
    colors = ["green", "red", "blue", "yellow", "magenta"]
    plots = []
    for name in scenario_names:
        results = all_cumulative_cases_per_day[name]
        for r in results:
            cumulative_cases = []
            for d in range(num_days):
                cumulative_cases.append(r[d])
            p, = plt.plot(dates, cumulative_cases, color=colors[i])
        plots.append(p)
        i += 1

    plt.xlabel("Day")
    plt.xticks(dates[::5], [x.strftime("%d/%m") for x in dates[::5]], rotation=45)
    plt.ylabel("Cumulative cases")
    plt.legend(plots, scenario_display_names)
    plt.tight_layout()
    plt.show()


    """
    plt.savefig("new_cases_per_day_comparison.png")
    plt.clf()

    i = 0
    for name in scenario_names:
        results = all_new_cases_per_day[name]
        days = range(num_days)
        means = []
        for d in days:
            results_for_day = []
            for r in results:
                results_for_day.append(r[d])
            means.append(sum(results_for_day) / len(results_for_day))
        plt.plot(days, means, color=colors[i])
        i += 1

    plt.xlabel("Day")
    plt.xticks(range(num_days)[::5], [x.strftime("%d/%m") for x in dates[::5]], rotation=45)
    plt.ylabel("Mean number of new infections")
    plt.legend(scenario_display_names)
    plt.tight_layout()
    plt.savefig("mean_new_cases_per_day_comparison.png")
    plt.clf()


    # FIXME
    #scenario_display_names[2] = "Geographically\n clustered non-compliers"

    #scenario_display_names = ["Basline", "25% \nnon-compliers", "50% \nnon-compliers", "75% \nnon-compliers", "100% \nnon-compliers"]
    #scenario_display_names = ["Baseline", "33K \nnon-compliers", "67K \nnon-compliers", "100K \nnon-compliers", "134K \nnon-compliers"]
    plt.boxplot(all_attack_rates, labels=scenario_display_names)
    plt.xticks(rotation=45)
    plt.ylabel("Attack rate over {} days".format(num_days))
    plt.tight_layout()

    plt.savefig("attack_rates_comparison.png")
    plt.clf()

    plt.boxplot(all_nums_non_compliers, labels=scenario_display_names)
    plt.xticks(rotation=45)
    plt.ylabel("Total number of non-compliers")
    plt.tight_layout()
    plt.savefig("num_non_compliers_comparison.png")
    plt.clf()"""


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", type=str)
    parser.add_argument("scenario_names", type=str, nargs="+")
    parser.add_argument("--scenario_display_names", type=str, nargs="+", default=None)
    args = parser.parse_args()
    main(args.output_dir, args.scenario_names, args.scenario_display_names)
