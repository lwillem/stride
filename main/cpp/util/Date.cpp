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
 * Implementation of the Date class.
 */

#include "Date.h"

#include <stdexcept>

namespace stride {
namespace util {

namespace {

// Howard Hinnant's date algorithms (public domain):
// http://howardhinnant.github.io/date_algorithms.html
// Converts a proleptic Gregorian (y, m, d) date to days since 1970-01-01.
long DaysFromCivil(int y, unsigned m, unsigned d)
{
        y -= m <= 2 ? 1 : 0;
        const long     era = (y >= 0 ? y : y - 399) / 400;
        const unsigned yoe = static_cast<unsigned>(y - era * 400);                    // [0, 399]
        const unsigned doy = (153 * (m + (m > 2 ? -3 : 9)) + 2) / 5 + d - 1;           // [0, 365]
        const unsigned doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;                    // [0, 146096]
        return era * 146097 + static_cast<long>(doe) - 719468;
}

// Inverse of DaysFromCivil: days since 1970-01-01 to proleptic Gregorian (y, m, d).
void CivilFromDays(long z, int& y, unsigned& m, unsigned& d)
{
        z += 719468;
        const long     era = (z >= 0 ? z : z - 146096) / 146097;
        const unsigned doe = static_cast<unsigned>(z - era * 146097);                 // [0, 146096]
        const unsigned yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;    // [0, 399]
        const long     yr  = static_cast<long>(yoe) + era * 400;
        const unsigned doy = doe - (365 * yoe + yoe / 4 - yoe / 100);                  // [0, 365]
        const unsigned mp  = (5 * doy + 2) / 153;                                      // [0, 11]
        d                  = doy - (153 * mp + 2) / 5 + 1;                             // [1, 31]
        m                  = mp + (mp < 10 ? 3 : static_cast<unsigned>(-9));           // [1, 12]
        y                  = static_cast<int>(yr) + (m <= 2 ? 1 : 0);
}

} // namespace

Date::Date(int year, unsigned month, unsigned day) : m_days(DaysFromCivil(year, month, day)) {}

Date Date::FromString(const std::string& s)
{
        if (s.size() != 10 || s[4] != '-' || s[7] != '-') {
                throw std::runtime_error("stride::util::Date::FromString> invalid date '" + s + "', expected YYYY-MM-DD");
        }
        const int      year  = std::stoi(s.substr(0, 4));
        const unsigned month = static_cast<unsigned>(std::stoi(s.substr(5, 2)));
        const unsigned day   = static_cast<unsigned>(std::stoi(s.substr(8, 2)));
        return Date(year, month, day);
}

Date Date::operator+(int days) const { return Date(m_days + days); }

int Date::operator-(const Date& other) const { return static_cast<int>(m_days - other.m_days); }

int Date::Year() const
{
        int      y;
        unsigned m, d;
        CivilFromDays(m_days, y, m, d);
        return y;
}

unsigned Date::Month() const
{
        int      y;
        unsigned m, d;
        CivilFromDays(m_days, y, m, d);
        return m;
}

unsigned Date::Day() const
{
        int      y;
        unsigned m, d;
        CivilFromDays(m_days, y, m, d);
        return d;
}

unsigned Date::DayOfWeek() const
{
        return static_cast<unsigned>(m_days >= -4 ? (m_days + 4) % 7 : (m_days + 5) % 7 + 6);
}

} // namespace util
} // namespace stride
