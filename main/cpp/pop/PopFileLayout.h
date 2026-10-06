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
 *  Copyright 2026, Willem L
 */

/**
 * @file
 * Column layout of a population file.
 */

#pragma once

#include <array>
#include <string>
#include <vector>

namespace stride {

/// One person as read from a row of the population file.
struct PopRecord
{
        unsigned int person_id          = 0U;
        unsigned int age                = 0U;
        unsigned int profession         = 0U;
        unsigned int household          = 0U;
        unsigned int school             = 0U;
        unsigned int workplace          = 0U;
        unsigned int community_weekend  = 0U;
        unsigned int community_weekday  = 0U;
        unsigned int household_cluster  = 0U;
        unsigned int collectivity       = 0U;
};

/// Which column of a population file holds which field. Built once from the first line of
/// the file, then used to parse every data row. Shared by PopBuilder (reading persons) and
/// PopSnapshotWriter (mirroring the input layout in its output).
class PopFileLayout
{
public:
        enum class Field
        {
                Age,
                PersonId,
                Profession,
                Household,
                School,
                Workplace,
                CommunityWeekend,
                CommunityWeekday,
                HouseholdCluster,
                Collectivity,
                Count
        };

        /// Determine the layout from the first line of a population file.
        static PopFileLayout FromFirstLine(const std::string& firstLine);

        /// Column separator: ";" if the first line contains one, "," otherwise.
        const std::string& Separator() const { return m_separator; }

        /// Whether the first line is a header (and not the first data row).
        bool HasHeader() const { return m_has_header; }

        /// Whether the file has a column for this field.
        bool Has(Field f) const { return m_column[static_cast<std::size_t>(f)] >= 0; }

        /// Parse one data row. Fields without a column keep their default (0), except the
        /// person id, which defaults to defaultPersonId.
        PopRecord Parse(const std::string& line, unsigned int defaultPersonId) const;

private:
        PopFileLayout();

        void Set(Field f, int column) { m_column[static_cast<std::size_t>(f)] = column; }

        std::string                                             m_separator;
        bool                                                    m_has_header;
        std::size_t                                             m_num_columns; ///< in the first line
        std::array<int, static_cast<std::size_t>(Field::Count)> m_column;      ///< -1 if absent
};

} // namespace stride
