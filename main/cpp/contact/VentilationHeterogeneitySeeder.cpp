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
 *  Copyright 2023
 */

/**
 * @file
 * Implementation for the Ventilation Heterogeneity class and Ventilation non compliance.
 */

#include "VentilationHeterogeneitySeeder.h"

#include "contact/ContactType.h"
#include "pop/Population.h"
#include "contact/ContactPoolSys.h"
#include "util/Exception.h"
#include "util/FileSys.h"
#include "util/RnMan.h"

#include <boost/property_tree/ptree.hpp>
#include <cassert>

using namespace boost::property_tree;
using namespace stride::util;
using namespace std;

namespace stride {

using namespace ContactType;

VentilationHeterogeneitySeeder::VentilationHeterogeneitySeeder(const ptree& config, RnMan& rnMan) : m_config(config), m_rn_man(rnMan) {}

shared_ptr<Population> VentilationHeterogeneitySeeder::Ventilation(shared_ptr<Population> pop)
{
	auto& population = *pop;

	auto& logger = population.RefEventLogger();

	boost::optional<string> ventilationFile = m_config.get_optional<string>("run.ventilation_file");    
    if (ventilationFile) {

        const auto fp = m_config.get<bool>("run.use_install_dirs") ? FileSys::GetDataDir() / (ventilationFile.value_or("")) : filesys::path(ventilationFile.value_or(""));
        ptree ventilationPt = FileSys::ReadPtreeFile(fp);
		
		for (ContactType::Id typ : ContactType::IdList) {
			std::string typString = ToString(typ);
			double ventilationInfo = ventilationPt.get<double>("ventilation_reduction." + typString,0);
			for (auto& pool: pop->RefPoolSys().RefPools(typ)) {
				pool.SetVentilation(ventilationInfo);
				logger->info("[VEN] {} {} {}", typString, ventilationInfo, pool.GetVentilation());
			};
		};


        } 
        
	return pop;

}

shared_ptr<Population> VentilationHeterogeneitySeeder::Seed(shared_ptr<Population> pop)
{

	auto& population = *pop;

	auto& logger = population.RefEventLogger();

	// Seed non-compliance

	// Non-compliance in pools
	boost::optional<string> nonCompliancePooltype = m_config.get_optional<string>("run.non_compliance_pooltype");
	if (nonCompliancePooltype) {
		string nonComplianceType = m_config.get<string>("run.non_compliance_type");
		ContactType::Id nonComplianceTypeId = ToId(nonComplianceType);
		for (auto& pool:pop->RefPoolSys().RefPools(nonComplianceTypeId)){
			pool.SetNonComplier();
			logger->info("[CHG] {} ", pool.IsNonComplier());
		}
	}

	return pop;
}



} // namespace stride
