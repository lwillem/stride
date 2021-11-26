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
    Some utility functions for running the superspreading scripts.
"""

import csv
import multiprocessing
import numpy as np
import os
import random
import subprocess
import xml.etree.ElementTree as ET

from collections import Counter
from scipy.stats import gamma
from scipy.integrate import quad

from postprocessing_util import get_output


def create_config(scenario_name, exp_id, contact_distribution, contact_distribution_overdispersion,
                    disease_config_file, event_log_level, holidays_file,
                    num_days, num_infected_seeds, population_file,
                    rng_seed, run_simplified, start_date, track_index_case,
                    transmission_probability_distribution, transmission_probability,
                    transmission_probability_distribution_overdispersion,
                    infected_seed_id=None):
    config = {
        "exp_id": exp_id,

        "age_contact_matrix_file": "contact_matrix_flanders_conditional_teachers.xml",
        "case_finding_capacity": 0,
        "cnt_intensity_householdCluster": 0,
        "contact_distribution": contact_distribution,
        "contact_distribution_overdispersion": contact_distribution_overdispersion,
        "delay_isolation_index": 0,
        "delay_contact_tracing": 0,
        "detection_probability": 0,
        "disease_config_file": disease_config_file,
        "event_log_level": event_log_level,
        "event_output_file": "true",
        "holidays_file": holidays_file,
        "hosp_probability_factor": 1,
        "immunity_profile": "None",
        "immunity_rate": 0,
        "num_cea_samples": 10000,
        "num_daily_imported_cases": 0,
        "num_days": num_days,
        "num_infected_seeds": num_infected_seeds,
        "num_participants_survey": 0,
        "num_threads": 1,
        "output_cases": "false",
        "output_persons": "false",
        "output_prefix": os.path.join("sim_output", scenario_name, "exp{:04}".format(exp_id)),
        "output_summary": "true",
        "population_file": population_file,
        "population_type": "default",
        "rng_seed": rng_seed,
        "run_simplified": run_simplified,
        "run_tag": scenario_name,
        "r0": 0,
        "seeding_age_max": 99,
        "seeding_age_min": 1,
        "start_date": start_date,
        "stride_log_level": "info",
        "test_false_negative": 0,
        "tracing_efficiency_household": 0,
        "tracing_efficiency_other": 0,
        "track_index_case": track_index_case,
        "transmission_probability_distribution": transmission_probability_distribution,
        "transmission_probability": transmission_probability,
        "transmission_probability_distribution_overdispersion": transmission_probability_distribution_overdispersion,
        "use_install_dirs": "true",
        "vaccine_link_probability": 0,
        "vaccine_profile": "None",
        "vaccine_rate": 0
    }

    if infected_seed_id is not None:
        config["infected_seed_id"] = infected_seed_id

    root = ET.Element("run")
    for name, value in config.items():
        parameter = ET.SubElement(root, name)
        parameter.text = str(value)

    tree = ET.ElementTree(root)
    tree.write("config/exp{:04}.xml".format(exp_id))

def f(t, shape, scale):
    return t * gamma.pdf(t, a=shape, scale=scale)

def get_mean_non_truncated_gamma(target_mean, shape):
    tolerance = 1.49e-4
    scale_est = target_mean / shape
    scale_params = np.arange(scale_est / 2, scale_est * 2, 0.0001)

    best_scale = np.nan
    best_mean = np.Inf

    for scale in scale_params:
        cdf1 = gamma.cdf(0, a=shape, scale=scale)
        cdf2 = gamma.cdf(1, a=shape, scale=scale)

        mean_tr = (quad(lambda x: f(x, shape, scale), 0, 1) / (cdf2 - cdf1))[0]

        if (abs(mean_tr - target_mean) < tolerance):
            if (abs(mean_tr - target_mean) < abs(best_mean - target_mean)):
                best_mean = mean_tr
                best_scale = scale

    return (best_scale * shape)

def run_and_summarize(exp_id):
    # Run
    subprocess.run(["./bin/stride", "-c exp{:04}.xml".format(exp_id)], stdout=subprocess.PIPE)
    # Summarize
    config = ET.parse('config/exp{:04}.xml'.format(exp_id)).getroot()
    output_prefix = (config.find('output_prefix').text).split("/")
    output_dir = output_prefix[0]
    scenario_name = output_prefix[1]

    output = get_output(output_dir, scenario_name, exp_id)

    num_days = output["parameters"]["num_days"]

    with open(os.path.join(output_dir, scenario_name, "exp{:04}".format(exp_id), "output_per_day.csv"), "w") as csvfile:
        fieldnames = ["sim_day", "num_cases", "r_t"]
        writer = csv.DictWriter(csvfile, fieldnames)
        writer.writeheader()

        for day in range(num_days):
            num_cases = 0
            r_t = np.nan
            if day in output["cases_per_day"]:
                num_cases = output["cases_per_day"][day]
            if day in output["rt_by_day"]:
                r_t = output["rt_by_day"][day]
            writer.writerow({
                "sim_day": day,
                "num_cases": num_cases,
                "r_t": r_t,
            })

    with open(os.path.join(output_dir, scenario_name, "exp{:04}".format(exp_id), "secondary_cases_frequencies.csv"), "w") as csvfile:
        fieldnames = ["num_secondary_cases", "frequency"]
        writer = csv.DictWriter(csvfile, fieldnames)
        writer.writeheader()

        frequencies = Counter(output["secondary_cases_by_individual"].values())
        for num_secondary_cases, frequency in frequencies.items():
            writer.writerow({
                "num_secondary_cases": num_secondary_cases,
                "frequency": frequency
            })

    with open(os.path.join(output_dir, scenario_name, "exp{:04}".format(exp_id), "secondary_cases_by_index_case.csv"), "w") as csvfile:
        fieldnames = ["index_case_id", "num_secondary_cases"]
        writer = csv.DictWriter(csvfile, fieldnames)
        writer.writeheader()

        for index_case_id in output["index_case_ids"]:
            writer.writerow({
                "index_case_id": index_case_id,
                "num_secondary_cases": output["secondary_cases_by_individual"][index_case_id]
            })

    with open(os.path.join(output_dir, scenario_name, "exp{:04}".format(exp_id), "transmissions_by_location.csv"), "w") as csvfile:
        fieldnames = ["location", "num_cases"]
        writer = csv.DictWriter(csvfile, fieldnames)
        writer.writeheader()

        for location, num_cases in output["transmissions_by_location"].items():
            writer.writerow({
                "location": location,
                "num_cases": num_cases
            })

    with open(os.path.join(output_dir, scenario_name, "exp{:04}".format(exp_id), "summary.csv"), "r") as csvinput:
        with open(os.path.join(output_dir, scenario_name, "exp{:04}".format(exp_id), "output_summary.csv"), "w") as csvoutput:
            writer = csv.writer(csvoutput, lineterminator='\n')
            reader = csv.reader(csvinput)

            all_rows = []

            # header
            row = next(reader)
            row.append("P80")
            all_rows.append(row)

            for row in reader:
                row.append(output["p80"])
                all_rows.append(row)

            writer.writerows(all_rows)

    # Copy config file to output dir
    subprocess.run(["mv", "config/exp{:04}.xml".format(exp_id),
                    os.path.join(output_dir, scenario_name, "exp{:04}".format(exp_id), "config.xml")])

    # Delete files that are no longer needed
    subprocess.run(["rm", os.path.join(output_dir, scenario_name, "exp{:04}".format(exp_id), "event_log.txt")])
    subprocess.run(["rm", os.path.join(output_dir, scenario_name, "exp{:04}".format(exp_id), "summary.csv")])

def run_parallel(scenario_name, contact_distribution, contact_distribution_overdispersion,
                    disease_config_file, event_log_level, holidays_file, num_days,
                    num_infected_seeds, population_file, run_simplified, start_date,
                    track_index_case, transmission_probability_distribution,
                    transmission_probability, transmission_probability_distribution_overdispersion,
                    num_runs, infected_seed_id=None, num_parallel_workers=4):

    ########################################
    # Create config files for experiments. #
    ########################################

    if not os.path.isdir("sim_output"):
        os.mkdir("sim_output")
    if not os.path.isdir(os.path.join("sim_output", scenario_name)):
        os.mkdir(os.path.join("sim_output", scenario_name))

    rngs = [random.randrange(10**9) for i in range(num_runs)]
    with open(os.path.join("sim_output", scenario_name, "exp_ids.txt"), "w") as f:
        for exp_id in range(1, num_runs + 1):
            print(exp_id, file=f)
        #f.writelines([str(exp_id) for exp_id in range(1, num_runs + 1)])
    for i in range(1, num_runs + 1):
        create_config(scenario_name=scenario_name, exp_id=i,
                    contact_distribution=contact_distribution, contact_distribution_overdispersion=contact_distribution_overdispersion,
                    disease_config_file=disease_config_file, event_log_level=event_log_level,
                    holidays_file=holidays_file, num_days=num_days, num_infected_seeds=num_infected_seeds,
                    population_file=population_file, rng_seed=rngs[i-1],
                    run_simplified=run_simplified, start_date=start_date, track_index_case=track_index_case,
                    transmission_probability_distribution=transmission_probability_distribution,
                    transmission_probability=transmission_probability,
                    transmission_probability_distribution_overdispersion=transmission_probability_distribution_overdispersion,
                    infected_seed_id=infected_seed_id)

    ######################################################
    # Run simulations + some postprocessing in parallel. #
    ######################################################

    with multiprocessing.Pool(processes=num_parallel_workers) as pool:
        pool.map(run_and_summarize, list(range(1, num_runs + 1)))
