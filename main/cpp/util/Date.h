/*
 *  This is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *  The software is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  You should have received a copy of the GNU General Public License
 *  along with the software. If not, see <http://www.gnu.org/licenses/>.
 *
 *  Copyright 2026
 */

/**
 * @file
 * A minimal Gregorian calendar date, replacing boost::gregorian::date for
 * the subset stride needs: parsing "YYYY-MM-DD", adding/subtracting whole
 * days, and extracting year/month/day/day-of-week. No timezones, no locales.
 */

#pragma once

#include <string>

namespace stride {
namespace util {

class Date
{
public:
        Date() = default;
        Date(int year, unsigned month, unsigned day);

        /// Parse a date formatted as "YYYY-MM-DD".
        static Date FromString(const std::string& s);

        /// Add (or subtract, if negative) a number of days.
        Date operator+(int days) const;

        /// Number of days between two dates (this minus other).
        int operator-(const Date& other) const;

        bool operator==(const Date& other) const { return m_days == other.m_days; }
        bool operator!=(const Date& other) const { return m_days != other.m_days; }
        bool operator<(const Date& other) const { return m_days < other.m_days; }
        bool operator<=(const Date& other) const { return m_days <= other.m_days; }
        bool operator>(const Date& other) const { return m_days > other.m_days; }
        bool operator>=(const Date& other) const { return m_days >= other.m_days; }

        int Year() const;
        unsigned Month() const;
        unsigned Day() const;

        /// Day of week: 0 = Sunday, ..., 6 = Saturday (matches boost::gregorian's numbering).
        unsigned DayOfWeek() const;

private:
        explicit Date(long daysSinceEpoch) : m_days(daysSinceEpoch) {}

        long m_days = 0; ///< Days since 1970-01-01 (proleptic Gregorian calendar).
};

} // namespace util
} // namespace stride
