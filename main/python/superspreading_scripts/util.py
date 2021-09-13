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
#  Copyright 2020, Willem L, Kuylen E & Broeckhove J
############################################################################ #

"""
    Some utility functions for superspreading postprocessing scripts.
"""

import csv
import matplotlib.pyplot as plt
import numpy as np
import os

def get_experiment_ids(output_dir, scenario_name):
    exp_design_file = os.path.join(output_dir, scenario_name, "exp_design.csv")
    experiment_ids = []
    with open(exp_design_file) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            exp_id = int(row["exp_id"])
            experiment_ids.append(exp_id)
    return experiment_ids

def get_cases_output(output_dir, scenario_name, experiment_id):
    print("Getting output for exp " + str(experiment_id))

    log_file = os.path.join(output_dir, scenario_name, "exp" + "{:04}".format(experiment_id), "event_log.txt")

    potential_infectors = {}
    cases_by_day = {}

    with open(log_file) as f:
        for line in f:
            line = line.split(" ")
            tag = line[0]
            if tag == "[PRIM]": # Index case
                infected_id = int(float(line[1]))
                if infected_id not in potential_infectors:
                    potential_infectors[infected_id] = 0
            elif tag == "[TRAN]": # Transmission
                infector_id = int(float(line[2]))
                infected_id = int(float(line[1]))

                # Add infected to potential infectors
                if infected_id not in potential_infectors:
                    potential_infectors[infected_id] = 0
                # Add infector to infectors (if not yet done)
                # And add this infection to total secondary cases count
                if infector_id not in potential_infectors:
                    potential_infectors[infector_id] = 1
                else:
                    potential_infectors[infector_id] += 1

                sim_day = int(line[6])
                if sim_day in cases_by_day:
                    cases_by_day[sim_day] += 1
                else:
                    cases_by_day[sim_day] = 1

    cases_output = {
        "secondary_cases_by_individual": potential_infectors,
        "cases_by_day": cases_by_day
    }
    return cases_output


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

def get_total_cases(cases_by_day, num_days):
    total_cases = 0
    for day, cases in cases_by_day.items():
        if day < num_days:
            total_cases += cases

    return total_cases

'''
def get_params_by_experiment(output_dir, scenario_name):
    exp_design_file = os.path.join(output_dir, scenario_name, "exp_design.csv")
    params = {}
    with open(exp_design_file) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            exp_id = int(row["exp_id"])
            params[exp_id] = {
                "transmission_probability" : float(row["transmission_probability"]),
                "transmission_probability_distribution": row["transmission_probability_distribution"],
                "transmission_probability_distribution_overdispersion": float(row["transmission_probability_distribution_overdispersion"]),

                "community_contact_distribution": row["community_contact_distribution"],
                "community_contact_distribution_overdispersion": float(row["community_contact_distribution_overdispersion"]),

                "num_infected_seeds": int(row["num_infected_seeds"]),
                "num_days": int(row["num_days"]),
            }

    return params

def get_trans_prob_by_exp(output_dir, scenario_name):
    experiments = {}
    summary_file = os.path.join(output_dir, scenario_name, scenario_name + "_summary.csv")
    with open(summary_file) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            exp_id = int(row["exp_id"])
            transmission_probability = float(row["transmission_probability"])
            experiments[exp_id] = transmission_probability

    return experiments

def get_cumulative_cases(output_dir, scenario_name, experiment_id, num_days, include_index_cases=False):
    log_file = os.path.join(output_dir, scenario_name, "exp" + "{:04}".format(experiment_id), "event_log.txt")
    cumulative_cases = 0
    with open(log_file) as f:
        for line in f:
            line = line.split(" ")
            tag = line[0]
            if tag == "[PRIM]" and include_index_cases:
                sim_day = int(line[6])
                if sim_day < num_days:
                    cumulative_cases += 1
            elif tag == "[TRAN]":
                sim_day = int(line[6])
                if sim_day < num_days:
                    cumulative_cases += 1
    return cumulative_cases


'''
