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
    Script to create holidays_file for scenario with period of social distancing.
"""

import argparse
import csv
import numpy as np

from datetime import date, datetime, timedelta

def main(start_date, use_holidays):
    num_days_pre_lockdown = 30
    num_days_lockdown = 60
    num_days_post_lockdown = 510

    pre_pandemic_contacts = { # Based on data for Belgium from Socrates platform
        "community": 5.871927,
        "workplace": 6.617716,
    }

    lockdown_contacts = { # Based on data for Belgium from CoMix study
        "community": [0.6560469, 0.6664077],
        "workplace": [0.2950804, 0.4317267],
    }

    post_lockdown_contacts = { # Based on data for Belgium from CoMix study
        "community": [1.468919, 1.525332, 1.574778, 2.081238, 2.024913, 1.384216],
        "workplace": [1.507625, 1.893812, 1.963227, 2.043224, 1.726014, 0.8268804,]
    }

    general_holidays = {
        2020: { 1: [1], 4: [12, 13], 5: [1, 21, 31], 6: [1], 7: [21], 11: [1, 11], 12: [25] },
        2021: { 1: [1], 4: [5], 5: [1, 13, 24], 7: [21], 8: [15], 11: [1, 11], 12: [25]},
    }

    general_holidays_dates = []
    for year in general_holidays:
        for month in general_holidays[year]:
            for day in general_holidays[year][month]:
                general_holidays_dates.append(date(year, month, day))

    school_holidays = {
        2020: {
            1: [1, 2, 3, 4, 5], 2: [24, 25, 26, 27, 28, 29], 4: [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19],
            7: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31],
            8: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31],
            11: [2, 3, 4, 5, 6, 7, 8], 12: [21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31]
        },
        2021: {
            1: [1, 2, 3], 2: [15, 16, 17, 18, 19, 20, 21], 4: [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18],
            7: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31],
            8: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31],
            11: [1, 2, 3, 4, 5, 6, 7], 12: [27, 28, 29, 30, 31]
        }
    }

    school_holidays_dates = []
    for year in school_holidays:
        for month in school_holidays[year]:
            for day in school_holidays[year][month]:
                school_holidays_dates.append(date(year, month, day))

    filename = "calendar_sprspr_social_distancing_2020_2021.csv"
    with open(filename, "w") as csvfile:
        fieldnames = ["category", "date", "value", "type", "age", "age_char"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()

        if use_holidays:
            # Add general holidays
            for day in general_holidays_dates:
                writer.writerow({
                    "category": "general",
                    "date": day,
                    "value": 1,
                    "type": "boolean",
                    "age": "NA",
                    "age_char": "NA"
                })
            # Add school holidays
            for day in school_holidays_dates:
                for age in range(26):
                    writer.writerow({
                        "category": "schools_closed",
                        "date": day,
                        "value": 1,
                        "type": "double",
                        "age": age,
                        "age_char": age
                    })

        # Add community distancing / workplace distancing / school closures during lockdown
        start_lockdown = datetime.strptime(start_date, "%Y-%m-%d") + timedelta(num_days_pre_lockdown)
        lockdown_community_distancing = 1 - (np.mean(lockdown_contacts["community"]) / pre_pandemic_contacts["community"])
        lockdown_workplace_distancing = 1 - (np.mean(lockdown_contacts["workplace"]) / pre_pandemic_contacts["workplace"])

        for i in range(num_days_lockdown):
            day = start_lockdown + timedelta(i)
            # Community distancing
            writer.writerow({
                "category": "community_distancing",
                "date": day,
                "value": lockdown_community_distancing,
                "type": "double",
                "age": "NA",
                "age_char": "NA"
            })
            # Workplace distancing
            writer.writerow({
                "category": "workplace_distancing",
                "date": day,
                "value": lockdown_workplace_distancing,
                "type": "double",
                "age": "NA",
                "age_char": "NA"
            })
            # School closures
            for age in range(26):
                writer.writerow({
                    "category": "schools_closed",
                    "date": day,
                    "value": 1,
                    "type": "double",
                    "age": age,
                    "age_char": age
                })

        # Add community distancing / workplace distancing post-lockdown
        end_lockdown = start_lockdown + timedelta(num_days_lockdown)
        post_lockdown_community_distancing = 1 - (np.mean(post_lockdown_contacts["community"]) / pre_pandemic_contacts["community"])
        post_lockdown_workplace_distancing = 1 - (np.mean(post_lockdown_contacts["workplace"]) / pre_pandemic_contacts["workplace"])

        for i in range(num_days_post_lockdown):
            day = end_lockdown + timedelta(i)
            # Community distancing
            writer.writerow({
                "category": "community_distancing",
                "date": day,
                "value": post_lockdown_community_distancing,
                "type": "double",
                "age": "NA",
                "age_char": "NA"
            })
            # Workplace distancing
            writer.writerow({
                "category": "workplace_distancing",
                "date": day,
                "value": post_lockdown_workplace_distancing,
                "type": "double",
                "age": "NA",
                "age_char": "NA"
            })

if __name__=="__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--start_date", type=str, default="2020-02-17")
    parser.add_argument("--no_holidays", dest="use_holidays", action="store_false")

    args = parser.parse_args()

    main(args.start_date, args.use_holidays)
