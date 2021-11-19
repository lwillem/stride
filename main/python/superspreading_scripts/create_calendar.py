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
    Script to create holidays_file for scenario with period of social distancing.
"""

import argparse
import csv

from datetime import date, timedelta

def main(lockdown_community_distancing, lockdown_workplace_distancing, post_lockdown_community_distancing, post_lockdown_workplace_distancing):
    start_date = date(2020, 2, 12)
    num_days_pre_lockdown = 30
    num_days_lockdown = 60
    num_days_post_lockdown = 110

    general_holidays = {
        2020: {
            1: [1],
            4: [12, 13],
            5: [1, 21, 31],
            6: [1],
            7: [21],
            11: [1, 11],
            12: [25]
        }
    }

    general_holidays_dates = []
    for year in general_holidays:
        for month in general_holidays[year]:
            for day in general_holidays[year][month]:
                general_holidays_dates.append(date(year, month, day))

    school_holidays = {
        2020: {
            1: [1, 2, 3, 4, 5],
            2: [24, 25, 26, 27, 28, 29],
            4: [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19],
            7: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31],
            8: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31],
            11: [2, 3, 4, 5, 6, 7, 8],
            12: [21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31]
        }
    }

    school_holidays_dates = []
    for year in school_holidays:
        for month in school_holidays[year]:
            for day in school_holidays[year][month]:
                    school_holidays_dates.append(date(year, month, day))

    file_name = "calendar_social_distancing"
    file_name = file_name + "_comm_" + str(int(lockdown_community_distancing * 100)) + "_" + str(int(post_lockdown_community_distancing * 100))
    file_name = file_name + "_work_" + str(int(lockdown_workplace_distancing * 100)) + "_" + str(int(post_lockdown_workplace_distancing * 100))
    file_name += ".csv"

    with open(file_name, "w") as csvfile:
        fieldnames = ["category", "date", "value", "type", "age", "age_char"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

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
        for day in general_holidays_dates:
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
        start_lockdown = start_date + timedelta(num_days_pre_lockdown)
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
            # School closure
            for age in range(26):
                writer.writerow({
                    "category": "schools_closed",
                    "date": day,
                    "value": 1,
                    "type": "double",
                    "age": age,
                    "age_char": age
                })


        end_lockdown = start_lockdown + timedelta(num_days_lockdown)
        # Add community distancing / workplace distancing after lockdown
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
    parser.add_argument("--lockdown_community_distancing", type=float, default=0.85)
    parser.add_argument("--lockdown_workplace_distancing", type=float, default=0.85)
    parser.add_argument("--post_lockdown_community_distancing", type=float, default=0.65)
    parser.add_argument("--post_lockdown_workplace_distancing", type=float, default=0.65)

    args = parser.parse_args()
    main(args.lockdown_community_distancing, args.lockdown_workplace_distancing, args.post_lockdown_community_distancing, args.post_lockdown_workplace_distancing)
