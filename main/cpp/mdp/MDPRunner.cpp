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
 * Implementation for MDPRunner.
 */

#include "MDPRunner.h"
#include "pop/Population.h"
#include "sim/Sim.h"

using namespace boost::property_tree;
using namespace std;

namespace stride {

MDPRunner::MDPRunner(const ptree& configPt, shared_ptr<Sim> sim)
        : SimController(configPt)
{
    m_simulator = std::move(sim);
}


} // namespace stride
