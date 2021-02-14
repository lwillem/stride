import argparse
import csv
import os

def main(output_dir, scenario_name):
    # Get experiment IDs
    experiment_ids = []
    summary_file = os.path.join(output_dir, scenario_name, scenario_name + "_summary.csv")
    with open(summary_file) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            experiment_ids.append(int(row["exp_id"]))

    # Count number of non-compliers in each experiment
    all_nums_non_compliers = []
    for exp_id in experiment_ids:
        num_non_compliers = 0
        log_file = os.path.join(output_dir, scenario_name, "exp" + "{:04}".format(exp_id), "contact_log.txt")
        with open(log_file) as f:
            for line in f:
                line = line.split()
                tag = line[0]
                if tag == "[NCOM]":
                    num_non_compliers += 1
        all_nums_non_compliers.append(num_non_compliers)

    avg_num_non_compliers = sum(all_nums_non_compliers) / len(all_nums_non_compliers)
    print(avg_num_non_compliers)



#['[NCOM]', '1713287', '4', '347758']

'''
OrderedDict([
    ('population_file', 'pop_belgium3000k_c500_teachers_censushh.csv'),
    ('num_days', '120'),
    ('population_size', '3000000'),
    ('seeding_rate', '3e-04'),
    ('r0', '2.5'),
    ('transmission_probability', '0.0602784'),
    ('immunity_rate', '0'),
    ('num_threads', '1'),
    ('rng_seed', '6'),
    ('run_time', '112969'),
    ('total_time', '112969'), ('num_cases', '432835'), ('AR', '0.144278'),
    ('output_prefix', 'sim_output/_geoclustering/exp0020'),
    ('start_date', '2020-02-28'), ('age_contact_matrix_file', 'contact_matrix_flanders_conditional_teachers.xml'), ('num_participants_survey', '0'), ('disease_config_file', 'disease_covid19_age.xml'), ('contact_log_level', 'Transmissions'), ('contact_output_file', 'true'), ('global_information_policy', 'NoGlobalInformation'), ('holidays_file', 'calendar_belgium_2020_covid19_may_workplace.json'), ('immunity_profile', 'None'), ('output_cases', 'true'), ('output_persons', 'false'),
    ('output_summary', 'true'), ('population_type', 'default'), ('seeding_age_max', '99'), ('seeding_age_min', '1'), ('stride_log_level', 'info'), ('track_index_case', 'false'), ('use_install_dirs', 'true'), ('vaccine_link_probability', '0'), ('vaccine_profile', 'None'), ('vaccine_rate', '0'), ('run_tag', '_geoclustering'), ('num_cea_samples', '10000'), ('adaptive_symptomatic_behavior', 'true'), ('cnt_reduction_work', '0'), ('age_break_school_types', '18'), ('cnt_reduction_other', '0.8'), ('compliance_delay', '14'), ('num_daily_imported_cases', '0'),
    ('cnt_reduction_work_exit', '0'), ('cnt_reduction_other_exit', '0'), ('non_compliance_type', 'Hotspots'),
    ('pools_in_hotspots_file', 'data/pop_belgium3000k_c500_teachers_censushh_households_in_hotspots.xml'),
    ('exp_id', '20')])
'''

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", type=str, help="Directory containing log output of simulations")
    parser.add_argument("scenario_name", type=str, help="Scenario name")
    args = parser.parse_args()
    main(args.output_dir, args.scenario_name)
