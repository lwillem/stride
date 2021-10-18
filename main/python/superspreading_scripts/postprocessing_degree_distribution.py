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
"""

import argparse
import csv
import matplotlib.pyplot as plt
import multiprocessing
import networkx as nx
import numpy as np
import os

from plots import save_figure
from util import get_experiment_ids

def get_degree_distribution(output_dir, scenario_name, experiment_id):
    G = nx.Graph()

    # Add all individuals in population to graph as nodes
    summary_file = os.path.join(output_dir, scenario_name, "exp" + "{:04}".format(experiment_id), "summary.csv")
    population_size = 0
    with open(summary_file) as csvfile:
        reader = csv.DictReader(csvfile)
        # Get first line
        row = next(reader)
        population_size = int(row["population_size"])

    for person_id in range(population_size):
        G.add_node(person_id)

    # Add edges between individuals that have contact
    log_file = os.path.join(output_dir, scenario_name, "exp" + "{:04}".format(experiment_id), "event_log.txt")
    number_of_edges = 0
    with open(log_file) as f:
        for line in f:
            line = line.split(" ")
            tag = line[0]
            if tag == "[CONT]":
                p1_id = int(float(line[1]))
                p2_id = int(float(line[17]))
                # Log contact in graph
                G.add_edge(p1_id, p2_id)
                number_of_edges += 1

    print(number_of_edges)

    # Calculate degree frequencies
    degree_freqs = nx.degree_histogram(G)

    #avg_degree = np.mean([degree[1] for degree in nx.degree(G)])
    #print(avg_degree)
    #print(nx.number_of_edges(G))

    # Normalize to total number of nodes
    degree_freqs = [freq / population_size for freq in degree_freqs]

    return degree_freqs

def main(output_dir, scenario_names, display_scenario_names):
    num_days = 1

    for s_i in range(len(scenario_names)):
        scenario_name = scenario_names[s_i]
        print(scenario_name)

        experiment_ids = get_experiment_ids(output_dir, scenario_name)
        with multiprocessing.Pool(processes=4) as pool:
            degree_dist = pool.starmap(get_degree_distribution, [(output_dir, scenario_name, exp_id) for exp_id in experiment_ids])
            for results in degree_dist:
                plt.bar(range(len(results)), results)
                plt.xlabel("Degree")
                plt.xlim(-0.5,55)
                plt.ylabel("Frequency")
                plt.ylim(0, 0.20)
                save_figure(output_dir, "degree_dist_" + scenario_name, extension="png", dpi=100)

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", type=str, help="Directory containing simulation results")
    parser.add_argument("scenario_names", type=str, nargs="+", help="Names of scenarios to be postprocessed")
    parser.add_argument("--display_scenario_names", type=str, nargs="+", default=[], help="Names for scenarios to be displayed on plots")

    args = parser.parse_args()
    main(args.output_dir, args.scenario_names, args.display_scenario_names)
