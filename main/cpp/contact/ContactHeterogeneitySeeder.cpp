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
 *  Copyright 2020, Kuylen E, Willem L, Broeckhove J
 */

/**
 * @file
 * Implementation for the NonComplianceSeeder class.
 */

#include "ContactHeterogeneitySeeder.h"

#include "pop/Population.h"
#include "util/Exception.h"
#include "util/FileSys.h"
#include "util/RnMan.h"

#include "util/Ptree.h"
#include <cassert>

using namespace stride::ContactType;
using namespace stride::util;
using namespace std;

namespace stride {

ContactHeterogeneitySeeder::ContactHeterogeneitySeeder(const ptree& config, std::shared_ptr<util::RnMan> rnMan) : m_config(config), m_rn_man(rnMan) {}

shared_ptr<Population> ContactHeterogeneitySeeder::Seed(shared_ptr<Population> pop)
{

	auto& population = *pop;

	const EventLogMode::Id log_level   = EventLogMode::ToMode(m_config.get<string>("run.event_log_level", "None"));
	auto& logger = population.RefEventLogger();

	boost::optional<string> contact_distribution = m_config.get_optional<string>("run.contact_distribution");

	if (contact_distribution) {

		// Get target overdispersion
		double contact_distribution_overdispersion = m_config.get<double>("run.contact_distribution_overdispersion");

		if (*contact_distribution == "Gamma") {
			// Use distribution with mean 1 and overdispersion = contact_distribution_overdispersion
			double shape = contact_distribution_overdispersion;
			double scale = 1 / shape;
			auto gamma_generator = m_rn_man->at(0U).GetGammaGenerator(shape, scale);

			// Seed community contact factors
			for (size_t i = 0; i < population.size(); ++i) {
				auto individual_contact_factor = gamma_generator();
				population[i].SetIndividualContactFactor(individual_contact_factor);

				// Log person details
				if (log_level == EventLogMode::Id::All) {
					logger->info("[CNTH] {} {}",
							population[i].GetId(), individual_contact_factor);
				}
			}

		}
	}

	return pop;
}

} // namespace stride
