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
 *  Copyright 2021
 */

/**
 * @file
 * Implementation for MDPRunner.
 */

#include "MDPRunner.h"
#include "pop/Population.h"
#include "sim/Sim.h"

using namespace stride::sim_event;
using namespace boost::property_tree;
using namespace std;

namespace stride {

MDPRunner::MDPRunner(const ptree& configPt, shared_ptr<Sim> sim)
        : SimRunner(configPt, sim), m_clock("total_clock"), m_config(configPt), m_sim(std::move(sim))
{
    Notify(Id::SetupBegin);
    m_clock.Start();
    Notify(Id::SetupEnd);
}

void MDPRunner::Step()
{
    // Prelims.
    m_clock.Start();

    // We are AtStart: no steps have taken yet, so signal AtStart.
    if (m_sim->GetCalendar()->GetSimulationDay() == 0) {
        Notify(Id::AtStart);
    }

    // Execute and signal Stepped
    m_sim->TimeStep();
    Notify(Id::Stepped);

    m_clock.Stop();
}

void MDPRunner::End()
{
    m_clock.Stop();
    Notify(Id::Finished);
    m_clock.Reset();
}

} // namespace stride
