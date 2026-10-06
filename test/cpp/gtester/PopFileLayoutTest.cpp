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
 * Unit tests for reading the column layout of a population file.
 */

#include "pop/PopFileLayout.h"

#include <gtest/gtest.h>
#include <stdexcept>

using namespace std;
using namespace stride;
using Field = PopFileLayout::Field;

namespace {

/// The six required fields, in the order age, household, school, work, weekend, weekday.
void ExpectBase(const PopRecord& r, unsigned int age, unsigned int hh, unsigned int school, unsigned int work,
                unsigned int weekend, unsigned int weekday)
{
        EXPECT_EQ(r.age, age);
        EXPECT_EQ(r.household, hh);
        EXPECT_EQ(r.school, school);
        EXPECT_EQ(r.workplace, work);
        EXPECT_EQ(r.community_weekend, weekend);
        EXPECT_EQ(r.community_weekday, weekday);
}

} // namespace

// The common layout, as in pop_belgium600k_c500_teachers_censushh.csv.
TEST(PopFileLayout, QuotedHeaderPrimarySecondary)
{
        const auto l = PopFileLayout::FromFirstLine(
            R"("age","household_id","school_id","work_id","primary_community","secondary_community")");
        EXPECT_TRUE(l.HasHeader());
        EXPECT_EQ(l.Separator(), ",");
        EXPECT_FALSE(l.Has(Field::PersonId));
        EXPECT_TRUE(l.PositionalColumns().empty());
        EXPECT_TRUE(l.IgnoredColumns().empty());
        const auto r = l.Parse("28,8285,0,13799,1,2", 7U);
        ExpectBase(r, 28, 8285, 0, 13799, 1, 2);
        EXPECT_EQ(r.person_id, 7U);
}

// As in pop_belgium10k_c500_teachers_censushh.csv: semicolons, person_id and worker, CRLF.
TEST(PopFileLayout, SemicolonWithProfessionAndCrlf)
{
        const auto l = PopFileLayout::FromFirstLine(
            "age;person_id;worker;household_id;school_id;work_id;primary_community;secondary_community\r");
        EXPECT_EQ(l.Separator(), ";");
        EXPECT_EQ(l.Column(Field::PersonId), 1);
        EXPECT_EQ(l.Column(Field::Profession), 2);
        const auto r = l.Parse("55;7217;1;4284;0;478;20;21\r", 0U);
        ExpectBase(r, 55, 4284, 0, 478, 20, 21);
        EXPECT_EQ(r.person_id, 7217U);
        EXPECT_EQ(r.profession, 1U);
        EXPECT_EQ(l.Parse("71;108;NA;108;0;0;1;1\r", 0U).profession, 0U);
}

// Synonyms: workplace_id and community_weekend/weekday, case-insensitive.
TEST(PopFileLayout, Synonyms)
{
        for (const auto* header : {"age,household_id,school_id,workplace_id,community_weekend,community_weekday",
                                   "Age,Household_ID,School_ID,Work_ID,Community_Weekend_ID,Community_Weekday_ID"}) {
                const auto l = PopFileLayout::FromFirstLine(header);
                EXPECT_TRUE(l.PositionalColumns().empty()) << header;
                ExpectBase(l.Parse("1,2,3,4,5,6", 0U), 1, 2, 3, 4, 5, 6);
        }
}

// Column order is no longer load-bearing.
TEST(PopFileLayout, ReorderedColumns)
{
        const auto l = PopFileLayout::FromFirstLine(
            "community_weekday,household_id,age,community_weekend,work_id,school_id");
        ExpectBase(l.Parse("6,2,1,5,4,3", 0U), 1, 2, 3, 4, 5, 6);
}

// household_cluster_id and collectivity_id may both be present, in any position.
TEST(PopFileLayout, ClusterAndCollectivityTogether)
{
        const auto l = PopFileLayout::FromFirstLine(
            "collectivity_id,age,household_id,school_id,work_id,primary_community,secondary_community,"
            "household_cluster_id\r");
        const auto r = l.Parse("9,1,2,3,4,5,6,8\r", 0U);
        ExpectBase(r, 1, 2, 3, 4, 5, 6);
        EXPECT_EQ(r.collectivity, 9U);
        EXPECT_EQ(r.household_cluster, 8U);
}

// The extra column at the end of a CRLF file: its name carries the trailing CR.
TEST(PopFileLayout, CrlfExtraColumn)
{
        const auto l = PopFileLayout::FromFirstLine(
            "age,household_id,school_id,work_id,primary_community,secondary_community,collectivity_id\r");
        EXPECT_TRUE(l.Has(Field::Collectivity));
        EXPECT_EQ(l.Parse("1,2,3,4,5,6,7\r", 0U).collectivity, 7U);
}

// A column with an unknown name is ignored, wherever it is.
TEST(PopFileLayout, UnknownColumnIgnored)
{
        const auto l = PopFileLayout::FromFirstLine(
            "age,district,household_id,school_id,work_id,primary_community,secondary_community,x_coord");
        ASSERT_EQ(l.IgnoredColumns().size(), 2U);
        ExpectBase(l.Parse("1,99,2,3,4,5,6,98", 0U), 1, 2, 3, 4, 5, 6);
}

// An unknown name where a required field would be positionally: read positionally, as before.
TEST(PopFileLayout, UnknownNamePositionalFallback)
{
        const auto l = PopFileLayout::FromFirstLine("age,hh,school_id,work_id,community_a,community_b");
        ASSERT_EQ(l.PositionalColumns().size(), 3U);
        ExpectBase(l.Parse("1,2,3,4,5,6", 0U), 1, 2, 3, 4, 5, 6);
}

// As in pop_belgium100k_c500_teachers_censushh.csv: person_id without worker. Positionally
// this file was read shifted by one column; by name it is read correctly.
TEST(PopFileLayout, PersonIdWithoutWorker)
{
        const auto l = PopFileLayout::FromFirstLine(
            "\"age\",\"person_id\",\"household_id\",\"school_id\",\"work_id\",\"primary_community\","
            "\"secondary_community\"\r");
        EXPECT_FALSE(l.Has(Field::Profession));
        const auto r = l.Parse("0,45,46,0,0,1,2\r", 0U);
        ExpectBase(r, 0, 46, 0, 0, 1, 2);
        EXPECT_EQ(r.person_id, 45U);
}

// No header: the first line is a person, read positionally; NA counts as a value.
TEST(PopFileLayout, NoHeader)
{
        for (const auto* line : {"56,1,0,0,3432,7468", "24,35896326,NA,1226,11,9", "\"56\",1,0,0,3432,7468\r"}) {
                const auto l = PopFileLayout::FromFirstLine(line);
                EXPECT_FALSE(l.HasHeader()) << line;
                EXPECT_EQ(l.Column(Field::CommunityWeekday), 5) << line;
        }
        const auto l = PopFileLayout::FromFirstLine("56,1,0,0,3432,7468,12");
        EXPECT_EQ(l.IgnoredColumns().size(), 1U);
        ExpectBase(l.Parse("56,1,0,0,3432,7468,12", 3U), 56, 1, 0, 0, 3432, 7468);
}

// A row shorter than the header leaves the missing optional fields at 0.
TEST(PopFileLayout, ShortRowDefaults)
{
        const auto l = PopFileLayout::FromFirstLine(
            "age,household_id,school_id,work_id,primary_community,secondary_community,household_cluster_id");
        EXPECT_EQ(l.Parse("1,2,3,4,5,6", 0U).household_cluster, 0U);
}

// Ambiguous or incomplete headers are errors, not silent misreads.
TEST(PopFileLayout, Errors)
{
        EXPECT_THROW(PopFileLayout::FromFirstLine(
                         "age,household_id,school_id,work_id,primary_community,community_weekend,secondary_community"),
                     runtime_error);
        EXPECT_THROW(PopFileLayout::FromFirstLine("household_id,age,school_id,work_id,primary_community"),
                     runtime_error);
}
