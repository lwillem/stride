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
 *  Copyright 2020 Willem L, Kuylen E, Stijven S & Broeckhove J
 */

/**
 * @file
 * Implementation of scenario tests running in batch mode.
 */

#include "ScenarioData.h"
#include "execs/SimController.h"
#include "pop/Population.h"
#include "sim/Sim.h"
#include "sim/SimRunner.h"
#include "util/ConfigInfo.h"
#include "util/FileSys.h"

#include <boost/property_tree/ptree.hpp>
#include <cmath>
#include <gtest/gtest.h>
#include <iomanip>
#include <sstream>
#include <tuple>

using namespace std;
using namespace stride;
using namespace stride::util;
using namespace ::testing;
using boost::property_tree::ptree;

namespace Tests {

class ScenarioRuns : public ::testing::TestWithParam<const char*>
{
public:
        /// TestCase set up.
        static void SetUpTestCase() {}

        /// Tearing down TestCase
        static void TearDownTestCase() {}

protected:
        /// Destructor has to be virtual.
        ~ScenarioRuns() override = default;

        /// Set up for the test fixture
        void SetUp() override {}

        /// Tearing down the test fixture
        void TearDown() override {}
};

class RunsDefault : public ScenarioRuns
{
};


void RunTest(const string& testTag, tuple<ptree, unsigned int, double> d, unsigned int numThreads)
{
        // -----------------------------------------------------------------------------------------
        // Scenario configuration and target numbers.
        // -----------------------------------------------------------------------------------------
        auto       config = get<0>(d);
        const auto target = get<1>(d);
        const auto margin = get<2>(d);
        config.put("run.num_threads", numThreads);
        config.put("run.output_prefix", string("tests/gtester_")
                   + ::testing::UnitTest::GetInstance()->current_test_info()->name() + "_" + testTag);

        // -----------------------------------------------------------------------------------------
        // Actual simulator run.
        // -----------------------------------------------------------------------------------------
        SimController controller(config, "TestController");
        controller.Control();
        const auto sim = controller.GetSim();
        const auto pop = sim->GetPopulation();

        // -----------------------------------------------------------------------------------------
        // Check results against target number (|res - target| < target * margin).
        // -----------------------------------------------------------------------------------------
        const unsigned int res = sim->GetPopulation()->GetTotalInfected();
        EXPECT_NEAR(res, target, target * margin)
            << "Failure at scenario: " << testTag << " with number of threads: " << numThreads << endl;


}

TEST_P(RunsDefault, defaultpop_single)
{
        const string testTag = GetParam();
        RunTest(testTag, ScenarioData::Get(testTag), 1U);
}

#ifdef _OPENMP
TEST_P(RunsDefault, defaultpop_multi)
{
        const string testTag = GetParam();
        RunTest(testTag, ScenarioData::Get(testTag), ConfigInfo::NumberAvailableThreads());
}
#endif


namespace {

const char* tags_influenza[] = {"influenza_a", "influenza_b", "influenza_c"};

const char* tags_measles[] = {"measles_16", "measles_26"};

const char* tags_r0[] = {"r0_0", "r0_4", "r0_8", "r0_12", "r0_16"};

const char* tags_covid19[] = {"covid19_base", "covid19_all", "covid19_daily", "covid19_distancing",
		"covid19_age_15min", "covid19_householdclusters", "covid19_tracing","covid19_tracing_all",
		"covid19_transm","covid19_transm_gamma",
		"covid19_suscept","covid19_suscept_age","covid19_suscept_adapt",
		"covid19_fitting","covid19_fitting_adapt"};

} // namespace

INSTANTIATE_TEST_CASE_P(covid19, RunsDefault, ValuesIn(tags_covid19));

INSTANTIATE_TEST_CASE_P(influenza, RunsDefault, ValuesIn(tags_influenza));

//INSTANTIATE_TEST_CASE_P(measles, RunsDefault, ValuesIn(tags_measles));

//INSTANTIATE_TEST_CASE_P(r0, RunsDefault, ValuesIn(tags_r0));


} // namespace Tests
