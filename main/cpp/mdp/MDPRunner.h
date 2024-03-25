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
 *  Copyright 2024
 */

/**
 * @file
 * Header for the MDPRunner class.
 */

#pragma once

#include "execs/SimController.h"

#include <boost/property_tree/ptree.hpp>
#include <memory>
#include <string>

namespace stride {

class Sim;
class Population;

/**
 * Based on SimController class.
 * The simulation runner drive simulator through time steps.
 * It's functions are:
 * \li invokes the simulator builder (@see SimulatorBuilder)
 * \li manages elapsed time clock
 * \li manages time steps
 */
class MDPRunner : public SimController
{
public:
        /// Initialization with property tree.
        /// \param configPt config info for run and for config of simulator
        explicit MDPRunner(const boost::property_tree::ptree& configPt, std::shared_ptr<Sim> sim);

};

} // namespace stride
