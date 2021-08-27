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
 *  Copyright 2020, Libin P, Willem L
 */

/**
 * @file
 * Implementation for the UniversalTesting class.
 */

#include "UniversalTesting.h"

#include "PublicHealthAgency.h"

#include "calendar/Calendar.h"
#include "pop/Population.h"
#include "util/Containers.h"
#include "util/CSV.h"
#include "util/FileSys.h"
#include "util/LogUtils.h"
#include "util/StringUtils.h"

#include <boost/property_tree/ptree.hpp>
#include <algorithm>
#include <stdexcept>

namespace stride {

using namespace boost::property_tree;
using namespace stride::ContactType;
using namespace stride::util;
using namespace std;

// Default constructor
UniversalTesting::UniversalTesting(): 
                m_unitest_planning_output_fn(),
	            m_unitest_pool_allocation(),
		        m_unitest_fnr(-1.0), m_unitest_n_tests_per_day(0), m_unitest_pool_size(0),
	            m_unitest_test_compliance(0.0), m_unitest_isolation_compliance(0.0),
                m_unitest_isolation_delay(0),
                m_unitest_detectable_delay(0),
                m_unitest_isolation_strategy("isolate-pool"),
	            m_unitest_planning(),
	            m_unitest_day_in_sweep(0)
	{}

void UniversalTesting::Initialize(const ptree& config){
	m_unitest_pool_allocation      = config.get<std::string>("run.unitest_pool_allocation","");
    m_unitest_fnr                  = config.get<double>("run.unitest_fnr",-1.0);
    m_unitest_n_tests_per_day      = config.get<unsigned int>("run.unitest_n_tests_per_day",0);
    m_unitest_pool_size            = config.get<unsigned int>("run.unitest_pool_size",0);
    m_unitest_test_compliance      = config.get<double>("run.unitest_test_compliance",0.0);
    m_unitest_isolation_compliance = config.get<double>("run.unitest_isolation_compliance",0.0);
    m_unitest_isolation_delay      = config.get<int>("run.unitest_isolation_delay",1);
    m_unitest_detectable_delay     = config.get<int>("run.unitest_detectable_delay",2);
    m_unitest_isolation_strategy   = config.get<std::string>("run.unitest_isolation_strategy","isolate-pool");
    const auto prefix = config.get<string>("run.output_prefix");
    m_unitest_planning_output_fn   = FileSys::BuildPath(prefix, "unitest_planning.csv");

    Replace(m_unitest_pool_allocation, "$unitest_pool_size", std::to_string(m_unitest_pool_size));
}

void UniversalTesting::PerformUniversalTesting(std::shared_ptr<Population> pop, 
        RnHandler& rnHandler, const std::shared_ptr<Calendar> calendar,
        PublicHealthAgency& pha)
{
  if (!calendar->IsUniversalTestingActivated())
    return;

  if (m_unitest_fnr < 0.0 || m_unitest_n_tests_per_day <= 0)
    return;

  const auto simDay = calendar->GetSimulationDay();

  if (m_unitest_planning.empty()) {
    const auto& households = pop->CRefPoolSys().CRefPools(Id::Household);

    std::map<std::string, std::map<int, PCRPool>> pools_per_georegion;
    unsigned int total_pools = 0;
    CSV allocation(m_unitest_pool_allocation);
    size_t georegion_idx = allocation.GetIndexForLabel("province");
    size_t pool_id_idx = allocation.GetIndexForLabel("pool_id");
    size_t household_id_idx = allocation.GetIndexForLabel("household_id");
    for (const auto& row : allocation) {
        std::string georegion = row.GetValue(georegion_idx);
        int pool_id = row.GetValue<int>(pool_id_idx);

        //if the georegion is not yet in the map, introduce it 
        if (pools_per_georegion.find(georegion) == pools_per_georegion.end()) {
            pools_per_georegion[georegion] = std::map<int,PCRPool>();
        }
        auto& pools = pools_per_georegion[georegion];

        //if the PCR pool is not yet in the map, introduce it 
        if (pools.find(pool_id) == pools.end()) {
            pools[pool_id] = PCRPool();
            pools[pool_id].SetId(pool_id);
            pools[pool_id].SetGeoRegion(georegion);

            ++total_pools;
        }
        auto& pool = pools[pool_id];

        int household_id = row.GetValue<int>(household_id_idx);
		const auto& household = households[household_id].GetPool();
        pool.AddHousehold(household);
    } 
    
    unsigned int n_days = ceil(total_pools / (float)m_unitest_n_tests_per_day);
    for (unsigned int day = 0; day < n_days; ++day) {
        m_unitest_planning.push_back(std::set<PCRPool>());
    }

    unsigned int day = 0;
    //TODO: shuffle regions keys randomly
    for (const auto& key_val: pools_per_georegion) {
        const auto& region = key_val.first;
        const auto& pools = util::MapValuesToVec(pools_per_georegion[region]);
        //TODO: pools can be shuffled randomly
        for (const auto& pool: pools) {
            m_unitest_planning[day].insert(pool);
            ++day;
            //when at the end of the planning, reset day
            if (day == m_unitest_planning.size())
               day = 0; 
        }
    }

    //write the planning to file
    ofstream of;
    of.open(m_unitest_planning_output_fn.c_str());
    of << "day,georegion,id,size" << std::endl;
    for (unsigned int day = 0; day < n_days; ++day) {
        for (const auto& pool: m_unitest_planning[day]) {
            of << day << "," 
               << pool.GetGeoRegion() << "," 
               << pool.GetId() << ","
               << pool.GetIndividuals().size()
               << std::endl;
        }
    }
    of.close();

#ifndef NDEBUG
    unsigned int enlisted_pop = 0;
    //test: check that the stride population and the population allocated over pools coincide
    for (unsigned int day = 0; day < n_days; ++day) {
        for (const auto& pool : m_unitest_planning[day]) {
            enlisted_pop += pool.GetIndividuals().size();
        }
    }
    std::cerr << "[UNIVERSAL] Enlisted vs total pop:"
        << enlisted_pop << " vs " << pop->size() << std::endl;
    assert(enlisted_pop == pop->size());

    //test: check that the pools' size does not exceed m_unitest_pool_size
    //test: report the nr of pools that are not completely full 
    //        (i.e., leftover_pools, there should not be many of them)
    int leftover_pools = 0;
    int filled_pools = 0;
    for (unsigned int day = 0; day < n_days; ++day) {
        for (const auto& pool : m_unitest_planning[day]) {
            assert(pool.GetIndividuals().size() <= m_unitest_pool_size);
            if (pool.GetIndividuals().size() < m_unitest_pool_size) {
                leftover_pools += 1;
            } else {
                filled_pools += 1;
            }
        }
    }
    std::cerr << "[UNIVERSAL] Leftover pools: " << leftover_pools << std::endl;
    std::cerr << "[UNIVERSAL] Filled pools: " << filled_pools << std::endl;

    //test: verify that the number of daily tests is not exceeded
    for (unsigned int day = 0; day < n_days; ++day) {
        std::cerr << "[UNIVERSAL]"
            << " Daily tests " << m_unitest_planning[day].size()
            << " on day " << day
            << " vs budget " <<  m_unitest_n_tests_per_day << std::endl;
        assert(m_unitest_planning[day].size() <= m_unitest_n_tests_per_day);
    }
#endif 
  }

  auto& logger = pop->RefEventLogger();

  //perform the testing, according to the planning
  for (const auto& pool : m_unitest_planning[m_unitest_day_in_sweep]) {
    logger->info("[UNITEST] {} {}", simDay, m_unitest_day_in_sweep);

    bool pool_positive_and_detectable = false;
    std::vector<std::vector<Person*>> tested_households;
    for (const auto& household : pool.GetHouseholds()) {
      bool compliant = rnHandler.Binomial(m_unitest_test_compliance);
     
      if (compliant) {
        tested_households.push_back(household);
      
        for (const Person* indiv : household) {
          auto h = indiv->GetHealth();
          if (h.IsInfected() && h.IsPcrDetectable(m_unitest_detectable_delay))
            pool_positive_and_detectable = true;
        }
      }
    }

    bool pcr_test_positive = pool_positive_and_detectable && rnHandler.Binomial(1-m_unitest_fnr);
    if (pcr_test_positive) {
      for (auto& household : tested_households) {
        if (m_unitest_isolation_strategy == "isolate-pool") {
          bool isolation_compliance = rnHandler.Binomial(m_unitest_isolation_compliance);
          if (isolation_compliance) {
            for (const auto& indiv : household) {
                unsigned int start = simDay + 1 + m_unitest_isolation_delay;
                indiv->Isolate(simDay, start, start + 7);
                logger->info("[UNITEST-ISOLATE] pool_id={} household_id={} indiv_id={} infected?={} isolation_delay={} sim_day={}",
                //logger->info("[UNITEST-ISOLATE] {} {} {} {} {} {}",
            	  	                               pool.GetId(),
                                                   indiv->GetPoolId(Id::Household),
                                                   indiv->GetId(), 
                                                   indiv->GetHealth().IsInfected(),
                                                   m_unitest_isolation_delay,
                                                   simDay);
            }
          }
        } else if (m_unitest_isolation_strategy == "trace") {
          for (Person* indiv : household) {
            auto h = indiv->GetHealth();
            if (h.IsInfected() && h.IsPcrDetectable(m_unitest_detectable_delay)) {
              bool pcr_test_positive = rnHandler.Binomial(1-m_unitest_fnr);
              if (pcr_test_positive)
                pha.Trace(*indiv, pop, rnHandler, calendar);
            }
          }
        } else {
          throw std::runtime_error("Invalid unitest_isolation_strategy");
        }
      }
    }
  }

  //move to the next day in the sweep,
  //if at the end of the sweep: reset
  ++m_unitest_day_in_sweep;
  if (m_unitest_day_in_sweep == m_unitest_planning.size()) {
    m_unitest_day_in_sweep = 0; 
  }
}

} // namespace stride
