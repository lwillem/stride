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
 * Header for the InfectorExec class.
 */

#pragma once

#include <spdlog/spdlog.h>

namespace stride {

class ContactPool;
class AgeContactProfile;
class TransmissionProfile;
class Population;

namespace util {

class RnHandler;

}

/// For use in the InfectorMap and Sim; executes infector.
typedef void(InfectorExec)(ContactPool& pool, const AgeContactProfile& profile,
                           const TransmissionProfile& trans_profile, util::RnHandler& rnHandler,
                           unsigned short int sim_day, std::shared_ptr<spdlog::logger> event_logger,
						   double cnt_reduction_work, double cnt_reduction_other, double cnt_reduction_school,
						   double collectivity_distancing_factor,
						   double cnt_reduction_intergeneration, unsigned int cnt_reduction_intergeneration_cutoff,
						   std::shared_ptr<Population> population, double m_cnt_intensity_householdCluster);

} // namespace stride
