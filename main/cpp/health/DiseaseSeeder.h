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
 *  Copyright 2018, Kuylen E, Willem L, Broeckhove J
 */

/**
 * @file
 * Header for the DiseaseSeeder class.
 */

#pragma once

#include "contact/ContactPool.h"
#include "TransmissionProfile.h"
#include "util/RnMan.h"
#include "util/SegmentedVector.h"

#include "util/Ptree.h"
#include <memory>

namespace stride {

class Population;
class Sim;

/**
 * Seeds population w.r.t infection.
 */
class DiseaseSeeder
{
public:
        /// Initializing DiseaseSeeder.
        DiseaseSeeder(const stride::util::ptree& config, std::shared_ptr<util::RnMan> rnMan);

        /// Import infected cases into the population
        void ImportInfectedCases(std::shared_ptr<Population> pop, unsigned int numInfected, unsigned int simDay, const TransmissionProfile& transProfile, util::Rn& rn);

private:
        const stride::util::ptree& m_config; ///< Run config.
        std::shared_ptr<util::RnMan>          m_rn_man; ///< Random number manager.
};

} // namespace stride
