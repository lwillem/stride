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
 * Implementation for SimRunner.
 */

#include "SimRunner.h"

#include "calendar/Calendar.h"
#include "pop/Population.h"
#include "sim/Sim.h"
#include "util/SummaryFile.h"

using namespace boost::property_tree;
using namespace std;

namespace stride {

SimRunner::SimRunner(const ptree& configPt, shared_ptr<Sim> sim)
    : m_clock("total_clock"), m_config(configPt), m_sim(std::move(sim))
{
        m_clock.Start();
}

void SimRunner::Run(unsigned int numSteps)
{
        if (numSteps != 0U) {

			m_clock.Start();
			const auto numDays = m_config.get<unsigned int>("run.num_days");

			// Take numSteps but do not go beyond numDays.
			for (unsigned int i = 0; i < numSteps && m_sim->GetCalendar()->GetSimulationDay() < numDays; i++) {
				   m_sim->TimeStep();
			}

			m_clock.Stop();
        }
}

void SimRunner::Run()
{
	Run(m_config.get<unsigned int>("run.num_days"));
}


void SimRunner::PrintSummary()
{
	const auto dur      = duration_cast<std::chrono::milliseconds>(GetClock().Get());
	const auto milli    = static_cast<unsigned int>(dur.count());

	SummaryFile  summary_file(m_config.get<string>("run.output_prefix"));

	summary_file.Print(m_config,
			static_cast<unsigned int>(m_sim->GetPopulation()->size()),
			m_sim->GetPopulation()->GetTotalInfected(),
			m_sim->RefTransmissionProfile().GetHomogeneousProbability(),
			milli, milli);
}

void SimRunner::Step()
{
    // Prelims.
    m_clock.Start();

    // Execute and signal Stepped
    m_sim->TimeStep();

    m_clock.Stop();
}

void SimRunner::End()
{
    m_clock.Stop();
    PrintSummary();
    m_clock.Reset();
}

} // namespace stride
