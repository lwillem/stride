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
 *  Copyright 2017, 2018, Kuylen E, Willem L, Broeckhove J
 */

/**
 * @file
 * Core Population class.
 */

#include "Population.h"

#include "disease/Health.h"
#include "util/Assert.h"
#include "util/FileSys.h"
#include "util/LogUtils.h"
#include "util/RnMan.h"
#include "util/RunConfigManager.h"
#include "util/SegmentedVector.h"
#include "util/StringUtils.h"

#include <boost/property_tree/ptree.hpp>
#include <utility>
#include "PopBuilder.h"

using namespace boost::property_tree;
using namespace std;
using namespace stride::util;
using namespace stride::ContactType;

namespace stride {

Population::Population() : m_pool_sys(), m_event_logger() {}

std::shared_ptr<Population> Population::Create(const boost::property_tree::ptree& config,
                                               std::shared_ptr<spdlog::logger> strideLogger)
{
        if (!strideLogger) {
                strideLogger = LogUtils::CreateNullLogger("Population_logger");
        }

        // --------------------------------------------------------------
        // Create empty population & and give it a InfectorLogger.
        // --------------------------------------------------------------
        const auto pop = Create();
        if (config.get<bool>("run.event_output_file", true)) {
                const auto prefix       = config.get<string>("run.output_prefix");
                const auto logPath      = FileSys::BuildPath(prefix, "event_log.txt");
                pop->RefEventLogger()   = LogUtils::CreateRotatingLogger("event_logger", logPath.string());
                pop->RefEventLogger()->set_pattern("%v");
                strideLogger->info("Event logging requested; logger set up.");
        } else {
                pop->RefEventLogger() = LogUtils::CreateNullLogger("event_logger");
                strideLogger->info("No Event logging requested.");
        }

        // -----------------------------------------------------------------------------------------
        // Build population.
        // -----------------------------------------------------------------------------------------
        strideLogger->info("Invoking PopBuilder.");
        PopBuilder(config, strideLogger).Build(pop);

        // -----------------------------------------------------------------------------------------
        // Done.
        // -----------------------------------------------------------------------------------------
        return pop;
}

std::shared_ptr<Population> Population::Create()
{

		// --------------------------------------------------------------
        // Create (empty) population and return it
        // --------------------------------------------------------------
        struct make_shared_enabler : public Population
        {
        };
        auto r = make_shared<make_shared_enabler>();
        return r;
}

Person* Population::CreatePerson(unsigned int id, double age, unsigned int householdId, unsigned int k12SchoolId,
                                 unsigned int college, unsigned int workId, unsigned int primaryCommunityId,
                                 unsigned int secondaryCommunityId, unsigned int householdClusterId, unsigned int collectivityId)
{
        return emplace_back(id, age, householdId, k12SchoolId, college, workId, primaryCommunityId,
                            secondaryCommunityId, householdClusterId, collectivityId);
}

unsigned int Population::GetTotalInfected() const
{
        unsigned int total{0U};
        for (const auto& p : *this) {
                const auto& h = p.GetHealth();
                total += h.IsInfected() || h.IsRecovered();
        }
        return total;
}

unsigned int Population::CountInfectedCases() const
{
        unsigned int total{0U};
        for (const auto& p : *this) {
                const auto& h = p.GetHealth();
                total += h.IsInfected();
        }
        return total;
}

unsigned int Population::CountExposedCases() const
{
        unsigned int total{0U};
        for (const auto& p : *this) {
                const auto& h = p.GetHealth();
                total += h.IsExposed();
        }
        return total;
}

unsigned int Population::CountInfectiousCases() const
{
        unsigned int total{0U};
        for (const auto& p : *this) {
                const auto& h = p.GetHealth();
                total += h.IsInfectious();
        }
        return total;
}

unsigned int Population::CountSymptomaticCases() const
{
        unsigned int total{0U};
        for (const auto& p : *this) {
                const auto& h = p.GetHealth();
                total += h.IsSymptomatic();
        }
        return total;
}

unsigned int Population::GetMaxAge() const
{
        unsigned int maxAge{0U};
        for (const auto& p : *this) {
                if(p.GetAge() > maxAge){
                	maxAge = p.GetAge();
                }
        }
        return maxAge;
}

unsigned int Population::GetPoolSize(ContactType::Id typeId, const Person* p) const {

	// get ContactPool id
	unsigned int poolId = p->GetPoolId(typeId);

	// poolId 0 means "not part of such a PoolType"
	if(poolId == 0){
		return 0;
	}

	// return ContactPool size
	return m_pool_sys.CRefPools(typeId)[poolId].size();
}


} // namespace stride
