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

from statsmodels.nonparametric.smoothers_lowess import lowess

def get_experiment_ids(output_dir, scenario_name):
    exp_design_file = os.path.join(output_dir, scenario_name, "exp_design.csv")
    experiment_ids = []
    with open(exp_design_file) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            exp_id = int(row["exp_id"])
            experiment_ids.append(exp_id)
    return experiment_ids

def get_parameters(output_dir, scenario_name, experiment_id):
    summary_file = os.path.join(output_dir, scenario_name, "exp" + "{:04}".format(experiment_id), "summary.csv")

    parameters = {}
    with open(summary_file) as csvfile:
        reader = csv.DictReader(csvfile)
        row = next(reader) # Get first line
        parameters["num_days"] = int(row["num_days"])
        parameters["population_size"] = int(row["population_size"])

    return parameters

def get_rt_by_day(infected_by_day, secondary_cases_by_individual, num_days):
    rt_by_day = {}

    for day in range(num_days):
        # Check if any individuals were infected on this day
        if day in infected_by_day and len(infected_by_day[day]) > 0:
            # Get number of secondary cases caused by each individual infected on this day
            secondary_cases_day = []
            for infector_id in infected_by_day[day]:
                secondary_cases_day.append(secondary_cases_by_individual[infector_id])
            rt_by_day[day] = np.mean(secondary_cases_day)
        else:
            rt_by_day[day] = np.nan

    return rt_by_day

def get_cases_output(output_dir, scenario_name, experiment_id):
    print("Getting output for exp " + str(experiment_id))

    # Get parameters from summary file
    parameters = get_parameters(output_dir, scenario_name, experiment_id)

    # Get output data from file
    log_file = os.path.join(output_dir, scenario_name, "exp" + "{:04}".format(experiment_id), "event_log.txt")

    potential_infectors = {}
    cases_by_day = {} # Keep track of number of cases per day
    infected_by_day = {} # Keep track of IDs of persons infected per day
    cnt_probabilities = []
    individual_contact_factors = []
    transmissions_by_location = {
        "Household": 0,
        "K12School": 0,
        "College": 0,
        "Workplace": 0,
        "PrimaryCommunity": 0,
        "SecondaryCommunity": 0,
    }
    contact_probabilities_by_location = {
        "Household": [],
        "K12School": [],
        "College": [],
        "Workplace": [],
        "PrimaryCommunity": [],
        "SecondaryCommunity": [],
    }

    total_transmissions = 0

    num_contacts_per_participant = {}

    with open(log_file) as f:
        for line in f:
            line = line.split(" ")
            tag = line[0]
            if tag == "[PRIM]": # Index case
                infected_id = int(float(line[1]))
                sim_day = int(line[6])
                if infected_id not in potential_infectors:
                    potential_infectors[infected_id] = 0

                if sim_day in cases_by_day:
                    cases_by_day[sim_day] += 1
                else:
                    cases_by_day[sim_day] = 1

                if sim_day in infected_by_day:
                    infected_by_day[sim_day].append(infected_id)
                else:
                    infected_by_day[sim_day] = [infected_id]

            elif tag == "[TRAN]": # Transmission
                total_transmissions += 1
                infector_id = int(float(line[2]))
                infected_id = int(float(line[1]))

                location = line[5]
                if location in transmissions_by_location:
                    transmissions_by_location[location] += 1
                else:
                    transmissions_by_location[location] = 1

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

                if sim_day in infected_by_day:
                    infected_by_day[sim_day].append(infected_id)
                else:
                    infected_by_day[sim_day] = [infected_id]
            elif tag == "[PART]":
                participant_id = int(float(line[1]))
                if participant_id not in num_contacts_per_participant:
                    num_contacts_per_participant[participant_id] = 0
            elif tag == "[CONT]":
                participant_id = int(float(line[1]))
                if participant_id in num_contacts_per_participant:
                    num_contacts_per_participant[participant_id] += 1
                else:
                    num_contacts_per_participant[participant_id] = 1
            elif tag == "[CNTH]":
                individual_contact_factors.append(float(line[2]))
            elif tag == "[CCNT]":
                cnt_probability = float(line[1])
                cnt_probabilities.append(cnt_probability)
                location = line[2].strip()
                if location in contact_probabilities_by_location:
                    contact_probabilities_by_location[location].append(cnt_probability)
                else:
                    contact_probabilities_by_location[location] = [cnt_probability]


    # Get P80
    # Get ...
    # Get Rt by day
    rt_by_day = get_rt_by_day(infected_by_day, potential_infectors, parameters["num_days"])

    cases_output = {
        "experiment_id": experiment_id,
        "parameters": parameters,
        "secondary_cases_by_individual": potential_infectors,
        "cases_by_day": cases_by_day,
        "infected_by_day": infected_by_day,
        "rt_by_day": rt_by_day,
        "contact_probabilities": cnt_probabilities,
        "individual_contact_factors": individual_contact_factors,
        "transmissions_by_location": transmissions_by_location,
        "contact_probabilities_by_location": contact_probabilities_by_location,
        "num_contacts_per_participant": num_contacts_per_participant
    }
    return cases_output

def get_day_of_last_infection(cases_by_day, num_days):
    day_of_last_infection = np.nan
    for day in range(num_days - 1, -1, -1):
        if day in cases_by_day and cases_by_day[day] > 0:
            day_of_last_infection = day
            break

    return day_of_last_infection

def get_num_cases_over_period(cases_by_day, start_day, end_day):
    total_cases = 0
    for day in range(start_day, end_day + 1):
        if day in cases_by_day:
            total_cases += cases_by_day[day]
    return total_cases

def get_herd_immunity_threshold(rt_by_day, cases_by_day, num_days, population_size, get_day=False):
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
            cumulative_cases += cases_by_day[day] if day in cases_by_day else 0
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

def get_total_cases(cases_by_day, num_days):
    total_cases = 0
    for day, cases in cases_by_day.items():
        if day < num_days:
            total_cases += cases

    return total_cases
