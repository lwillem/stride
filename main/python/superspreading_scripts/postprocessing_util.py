############################################################################ #
#  This file is part of the Stride software.
#  It is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or any
#  later version.
#  The software is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  You should have received a copy of the GNU General Public License,
#  along with the software. If not, see <http://www.gnu.org/licenses/>.
#  see http://www.gnu.org/licenses/.
#
#
#  Copyright 2021, Kuylen E
############################################################################ #

"""
    Some utility functions for superspreading postprocessing scripts.
"""

import csv
import networkx as nx
import numpy as np
import os

from statsmodels.nonparametric.smoothers_lowess import lowess

def get_experiment_ids(output_dir, scenario_name):
    exp_ids_file = os.path.join(output_dir, scenario_name, "exp_ids.txt")
    experiment_ids = []
    with open(exp_ids_file) as f:
        for line in f:
            experiment_ids.append(int(line))

    return experiment_ids

def get_day_of_last_infection(cases_per_day, num_days):
    day_of_last_infection = np.nan
    for day in range(num_days - 1, -1, -1):
        if day in cases_per_day and cases_per_day[day] > 0:
            day_of_last_infection = day
            break

    return day_of_last_infection

def get_degree_distribution(output_dir, scenario_name, experiment_id, population_size):
    G = nx.Graph()

    # Add all individuals in population to graph as nodes
    for person_id in range(population_size):
        G.add_node(person_id)

    # Add edges between individuals that have contact
    log_file = os.path.join(output_dir, scenario_name, "exp{:04}".format(experiment_id), "event_log.txt")
    with open(log_file) as f:
        for line in f:
            line = line.split(" ")
            tag = line[0]
            if tag == "[CONT]":
                p1_id = int(float(line[1]))
                p2_id = int(float(line[17]))
                # Add edge to graph
                G.add_edge(p1_id, p2_id)

    # Calculate degree frequencies
    degree_freqs = nx.degree_histogram(G)

    # Normalize to total number of nodes
    #degree_freqs = [freq / population_size for freq in degree_freqs]

    return degree_freqs

def get_degree_distribution_from_file(output_dir, scenario_name, experiment_id):
    output_file = os.path.join(output_dir, scenario_name, "exp{:04}".format(experiment_id), "degree_distribution.csv")
    degree_distribution = {}
    with open(output_file) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            degree = int(row["degree"])
            frequency = int(row["frequency"])
            degree_distribution[degree] = frequency

    return degree_distribution

def get_herd_immunity_threshold(rt_by_day, cases_per_day, num_days, population_size, get_day=False):
    total_cases = sum(cases_per_day.values())
    # Smooth (using LOWESS function)
    rt_by_day_smoothed = lowess([rt_by_day[day] for day in range(num_days)], range(num_days), is_sorted=True, return_sorted=False)

    # Look for last day where smoothed Rt >= 1
    last_day_rt_greq_1 = np.nan
    for day in range(num_days - 1, -1, -1):
        rt = rt_by_day_smoothed[day]
        if rt >= 1:
            last_day_rt_greq_1 = day
            break

    if get_day: # Return day on which herd immunity threshold is reached
        return last_day_rt_greq_1
    else:
        # Get proportion of population no longer susceptible on this day
        # = cumulative cases (including index cases) / population_size
        if np.isnan(last_day_rt_greq_1):
            return np.nan
        cumulative_cases = 0
        for day in range(last_day_rt_greq_1 + 1):
            cumulative_cases += cases_per_day[day] if day in cases_per_day else 0
        herd_immunity_threshold = cumulative_cases / population_size
        return herd_immunity_threshold

def get_p80(secondary_cases_by_individual, extinction_threshold = 0):
    total_cases = sum(list(secondary_cases_by_individual.values()))
    if total_cases >= extinction_threshold and total_cases > 0:
        secondary_cases_sorted = list(secondary_cases_by_individual.values())
        secondary_cases_sorted.sort(reverse=True)

        num_cases_responsible = 0
        num_cases_caused = 0
        for s in secondary_cases_sorted:
            num_cases_caused += s
            num_cases_responsible += 1
            if num_cases_caused >= (total_cases * 0.80):
                break

        p80 = num_cases_responsible / total_cases
        return p80
    else:
        return np.nan

def get_parameters(output_dir, scenario_name, experiment_id):
    summary_file = os.path.join(output_dir, scenario_name, "exp" + "{:04}".format(experiment_id), "summary.csv")

    parameters = {}
    with open(summary_file) as csvfile:
        reader = csv.DictReader(csvfile)
        row = next(reader) # Get first line
        parameters["num_days"] = int(row["num_days"])
        parameters["population_size"] = int(row["population_size"])
        parameters["transmission_probability"] = float(row["transmission_probability"])

    return parameters

def get_rt_by_day(ids_infected_by_day, secondary_cases_by_individual, num_days):
    rt_by_day = {}

    for day in range(num_days):
        # Check if any individuals were infected on this day
        if day in ids_infected_by_day and len(ids_infected_by_day[day]) > 0:
            # Get number of secondary cases caused by each individual infected on this day
            secondary_cases_day = []
            for infector_id in ids_infected_by_day[day]:
                secondary_cases_day.append(secondary_cases_by_individual[infector_id])
            rt_by_day[day] = np.mean(secondary_cases_day)
        else:
            rt_by_day[day] = np.nan

    return rt_by_day

def get_resurgence_probability(all_cases_per_day, start_lockdown, end_lockdown, num_days, resurgence_threshold):
    runs_with_cases_during_lockdown = 0
    runs_above_resurgence_threshold = 0

    for run in all_cases_per_day:
        num_cases_during_lockdown = get_cases_over_period(run, 30, 90)
        num_cases_after_lockdown = get_cases_over_period(run, 90, 600)
        if not num_cases_during_lockdown == 0:
            runs_with_cases_during_lockdown += 1
            if num_cases_after_lockdown >= resurgence_threshold:
                runs_above_resurgence_threshold += 1

    return runs_above_resurgence_threshold / runs_with_cases_during_lockdown

def get_total_cases(cases_per_day, num_days):
    total_cases = 0
    for day, cases in cases_per_day.items():
        if day < num_days:
            total_cases += cases

    return total_cases

def get_cases_over_period(cases_per_day, start_day, end_day):
    total_cases = 0
    for day in range(start_day, end_day):
        if day in cases_per_day:
            total_cases += cases_per_day[day]

    return total_cases

def get_output(output_dir, scenario_name, experiment_id):

    log_file = os.path.join(output_dir, scenario_name, "exp" + "{:04}".format(experiment_id), "event_log.txt")

    # Get parameters from summary file
    parameters = get_parameters(output_dir, scenario_name, experiment_id)

    cases_per_day = {}
    ids_infected_by_day = {}
    secondary_cases_by_individual = {}
    index_case_ids = []

    transmissions_by_location = {
        "Household": 0,
        "K12School": 0,
        "College": 0,
        "Workplace": 0,
        "PrimaryCommunity": 0,
        "SecondaryCommunity": 0,
    }

    with open(log_file) as f:
        for line in f:
            line = line.split(" ")
            tag = line[0]

            if tag == "[PRIM]": # Index case
                infected_id = int(float(line[1]))
                sim_day = int(line[6])
                if infected_id not in secondary_cases_by_individual:
                    secondary_cases_by_individual[infected_id] = 0

                index_case_ids.append(infected_id)

                if sim_day in cases_per_day:
                    cases_per_day[sim_day] += 1
                else:
                    cases_per_day[sim_day] = 1

                if sim_day in ids_infected_by_day:
                    ids_infected_by_day[sim_day].append(infected_id)
                else:
                    ids_infected_by_day[sim_day] = [infected_id]
            elif tag == "[TRAN]": # Transmission
                infector_id = int(float(line[2]))
                infected_id = int(float(line[1]))
                location = line[5]
                sim_day = int(line[6])

                # Add infected to potential infectors
                if infected_id not in secondary_cases_by_individual:
                    secondary_cases_by_individual[infected_id] = 0
                # Add infector to infectors (if not yet done)
                # And add this transmission to total secondary cases count
                if infector_id not in secondary_cases_by_individual:
                    secondary_cases_by_individual[infector_id] = 1
                else:
                    secondary_cases_by_individual[infector_id] += 1

                if sim_day in cases_per_day:
                    cases_per_day[sim_day] += 1
                else:
                    cases_per_day[sim_day] = 1

                if sim_day in ids_infected_by_day:
                    ids_infected_by_day[sim_day].append(infected_id)
                else:
                    ids_infected_by_day[sim_day] = [infected_id]

                if location in transmissions_by_location:
                    transmissions_by_location[location] += 1
                else:
                    transmissions_by_location[location] = 1

    rt_by_day = get_rt_by_day(ids_infected_by_day, secondary_cases_by_individual, parameters["num_days"])
    p80 = get_p80(secondary_cases_by_individual, extinction_threshold = 0)

    output = {
        "experiment_id": experiment_id,
        "parameters": parameters,
        "secondary_cases_by_individual": secondary_cases_by_individual,
        "index_case_ids": index_case_ids,
        "cases_per_day": cases_per_day,
        "rt_by_day": rt_by_day,
        "transmissions_by_location": transmissions_by_location,
        "p80": p80
    }

    return output

def get_summary_output(output_dir, scenario_name, exp_id, file_name="output_summary.csv"):
    summary_file = os.path.join(output_dir, scenario_name, "exp{:04}".format(exp_id), file_name)
    with open(summary_file) as csvfile:
        reader = csv.DictReader(csvfile)

        # Get first row
        row = next(reader)
        output = {
            "contact_matrix_file": row["age_contact_matrix_file"],
            "disease_config_file": row["disease_config_file"],
            "num_days": int(row["num_days"]),
            "population_size": int(row["population_size"]),
            "population_file": row["population_file"],
            "total_cases": int(row["num_cases"])
        }

        if "P80" in row:
            output["p80"] = float(row["P80"])

        return output

def get_num_secondary_cases_frequencies(output_dir, scenario_name, exp_id):
    secondary_cases_file = os.path.join(output_dir, scenario_name, "exp{:04}".format(exp_id), "secondary_cases_frequencies.csv")

    secondary_cases_frequencies = {}
    with open(secondary_cases_file) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            num_secondary_cases = int(row["num_secondary_cases"])
            frequency = int(row["frequency"])
            secondary_cases_frequencies[num_secondary_cases] = frequency

    return secondary_cases_frequencies

def get_num_secondary_cases_per_index_case(output_dir, scenario_name, exp_id):
    output_file = os.path.join(output_dir, scenario_name, "exp{:04}".format(exp_id), "secondary_cases_by_index_case.csv")
    all_secondary_cases = []
    with open(output_file) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            all_secondary_cases.append(int(row["num_secondary_cases"]))

    return np.mean(all_secondary_cases)

def get_output_per_day(output_dir, scenario_name, exp_id):
    num_cases_per_day = {}
    rt_by_day = {}
    output_file = os.path.join(output_dir, scenario_name, "exp{:04}".format(exp_id), "output_per_day.csv")
    with open(output_file) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            day = int(row["sim_day"])
            num_cases = int(row["num_cases"])
            r_t = float(row["r_t"])
            num_cases_per_day[day] = num_cases
            rt_by_day[day] = r_t

    return {"cases_per_day": num_cases_per_day, "rt_by_day": rt_by_day}

def get_transmissions_by_location(output_dir, scenario_name, exp_id):
    transmissions_by_location = {}

    output_file = os.path.join(output_dir, scenario_name, "exp{:04}".format(exp_id), "transmissions_by_location.csv")
    with open(output_file) as csvfile:
        reader = csv.DictReader(csvfile)

        for row in reader:
            location = row["location"]
            num_transmissions = int(row["num_cases"])
            transmissions_by_location[location] = num_transmissions

    return transmissions_by_location
