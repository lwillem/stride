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
///
/// - A first line without a header name (only numbers, empty or NA) is the first person;
///   the file is then read positionally: age, household, school, work, weekend, weekday.
/// - Otherwise columns are resolved by header name (case, quotes and a trailing CR are
///   ignored; synonyms such as work_id / workplace_id are accepted). A field without a
///   recognised name falls back to its positional column, provided that column's own name
///   was not recognised. Columns that remain unrecognised are ignored, not an error.
/// - An error is raised only for an ambiguous or incomplete header: two columns for one
///   field, or no column at all for a required field.
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
        bool Has(Field f) const { return Column(f) >= 0; }

        /// Zero-based column of this field, -1 if absent.
        int Column(Field f) const { return m_column[static_cast<std::size_t>(f)]; }

        /// Columns read by position because their name was not recognised (for logging).
        const std::vector<std::string>& PositionalColumns() const { return m_positional; }

        /// Columns that are not read at all (for logging).
        const std::vector<std::string>& IgnoredColumns() const { return m_ignored; }

        /// Canonical header name of a field.
        static std::string FieldName(Field f);

        /// Parse one data row. Fields without a column keep their default (0), except the
        /// person id, which defaults to defaultPersonId.
        PopRecord Parse(const std::string& line, unsigned int defaultPersonId) const;

private:
        PopFileLayout();

        /// The layout inferred from position only, as before columns were resolved by name.
        static PopFileLayout Positional(const std::vector<std::string>& names);

        /// Fields every population file must provide.
        static const std::vector<Field>& RequiredFields();

        void Set(Field f, int column) { m_column[static_cast<std::size_t>(f)] = column; }

        std::string                                             m_separator;
        bool                                                    m_has_header;
        std::array<int, static_cast<std::size_t>(Field::Count)> m_column; ///< -1 if absent
        std::vector<std::string>                                m_positional;
        std::vector<std::string>                                m_ignored;
};

} // namespace stride
