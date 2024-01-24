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
 * Header for the Simulation controller.
 */

#include "SimController.h"

#include "pop/Population.h"
#include "sim/Sim.h"
#include "util/ConfigInfo.h"
#include "util/FileSys.h"
#include "util/LogUtils.h"
#include "util/Stopwatch.h"
#include "util/SummaryFile.h"
#include "util/TimeStamp.h"

#include <boost/property_tree/xml_parser.hpp>
#include <regex>

using namespace std;
using namespace stride::util;
using namespace boost::property_tree;
using namespace boost::property_tree::xml_parser;

namespace stride {

SimController::SimController()
    : m_config(), m_output_prefix(), m_run_clock("run"),
	  m_stride_logger(nullptr), m_use_install_dirs(),
	  m_simulator(nullptr)

{
}

SimController::~SimController()
{
        Shutdown();
}

SimController::SimController(const ptree& config) : SimController()
{
        m_run_clock.Start();
        m_config           = config;
        m_output_prefix    = m_config.get<string>("run.output_prefix");
        m_use_install_dirs = m_config.get<bool>("run.use_install_dirs");
}

void SimController::CheckEnv()
{
        if (m_use_install_dirs) {
                auto log = [](const string& s) -> void { cerr << s << endl; };
                if (!FileSys::CheckInstallEnv(log)) {
                        throw runtime_error("SimController::CheckEnv> Install dirs not OK.");
                }
        }
}

void SimController::CheckOutputPrefix()
{
        if (FileSys::IsDirectoryString(m_output_prefix)) {
                FileSys::CreateDirectory(m_output_prefix);
        }
}

void SimController::InstallLogger()
{
        m_stride_logger = spdlog::get("stride_logger");
        // If there is no stride_logger yet ...
        if (!m_stride_logger) {
			const auto path = FileSys::BuildPath(m_output_prefix, "stride_log.txt");
			m_stride_logger = LogUtils::CreateCliLogger("stride_logger", path.string());
			spdlog::register_logger(m_stride_logger);
        }

        // spd log levels are: trace, debug, info, warn, error, critical, off
        const auto logLevel = m_config.get<string>("run.stride_log_level");
        m_stride_logger->set_level(spdlog::level::from_str(logLevel));
        m_stride_logger->flush_on(spdlog::level::err);

}

void SimController::Shutdown()
{
        m_run_clock.Stop();
        m_stride_logger->info("Shutting down after: {}", m_run_clock.ToString());
        m_stride_logger->flush();
        spdlog::drop("stride_logger");
}

void SimController::LogStartup()
{
        m_stride_logger->info("Starting up at: {}", TimeStamp().ToString());
        m_stride_logger->info("Executing revision: {}", ConfigInfo::GitRevision());
        m_stride_logger->info("Processor count: {}", ConfigInfo::ProcessorCount());
        m_stride_logger->info("Creating dir:  {}", m_output_prefix);
        m_stride_logger->trace("Executing:           {}", FileSys::GetExecPath().string());
        m_stride_logger->trace("Current directory:   {}", FileSys::GetCurrentDir().string());
        if (m_use_install_dirs) {
                m_stride_logger->trace("Install directory:   {}", FileSys::GetRootDir().string());
                m_stride_logger->trace("Config  directory:   {}", FileSys::GetConfigDir().string());
                m_stride_logger->trace("Data    directory:   {}", FileSys::GetDataDir().string());
        }
        if (ConfigInfo::HaveOpenMP()) {
                m_stride_logger->info("Max number OpenMP threads in this environment: {}",
                                      ConfigInfo::NumberAvailableThreads());
                m_stride_logger->info("Configured number of threads: {}",
                                      m_config.get<unsigned int>("run.num_threads"));
        } else {
                m_stride_logger->info("Not using OpenMP threads.");
        }
        stringstream ss;
        write_xml(ss, m_config, xml_writer_make_settings<ptree::key_type>(' ', 8));
        const auto s = ss.str();
        stringstream spretty;
        std::regex_replace(std::ostreambuf_iterator<char>(spretty), s.begin(), s.end(), std::regex("(\\n+)"), "\n");
        m_stride_logger->trace("Config :\n {}", spretty.str());
}

void SimController::Control()
{
        // -----------------------------------------------------------------------------------------
        // Prelims.
        // -----------------------------------------------------------------------------------------
        CheckEnv();
        CheckOutputPrefix();
        InstallLogger();
        LogStartup();

        // -----------------------------------------------------------------------------------------
        // Sim scenario: step 1, build a random number manager.
        // -----------------------------------------------------------------------------------------
        const RnInfo info{m_config.get<string>("run.rng_seed", "1,2,3,4"), "",
                          m_config.get<unsigned int>("run.num_threads")};
        RnMan        rnMan{info};

        // -----------------------------------------------------------------------------------------
        // Sim scenario: step 2, create a population, as described by the parameter in the config.
        // -----------------------------------------------------------------------------------------
        auto pop = Population::Create(m_config, m_stride_logger);

        // -----------------------------------------------------------------------------------------
        // Sim scenario: step 3, create a simulator, as described by the parameter in the config.
        // -----------------------------------------------------------------------------------------
        m_simulator = Sim::Create(m_config, pop, rnMan);

        // -----------------------------------------------------------------------------------------
        // Sim scenario: step 4, run and print results
        // -----------------------------------------------------------------------------------------
        Run();
        PrintSummary();

}

void SimController::Run(unsigned int numSteps)
{
        if (numSteps != 0U) {

        	m_run_clock.Start();
			const auto numDays = m_config.get<unsigned int>("run.num_days");

			// Take numSteps but do not go beyond numDays.
			for (unsigned int i = 0; i < numSteps && m_simulator->GetCalendar()->GetSimulationDay() < numDays; i++) {
				m_simulator->TimeStep();
			}

			m_run_clock.Stop();
        }
}

void SimController::Run()
{
	Run(m_config.get<unsigned int>("run.num_days"));
}


void SimController::PrintSummary()
{
	const auto  milli  = GetClock().ToUnsignedInteger();
	SummaryFile summary_file(m_config.get<string>("run.output_prefix"));

	summary_file.Print(m_config,
			static_cast<unsigned int>(m_simulator->GetPopulation()->size()),
			m_simulator->GetPopulation()->GetTotalInfected(),
			m_simulator->RefTransmissionProfile().GetHomogeneousProbability(),
			milli);
}

void SimController::Step()
{
    // Prelims.
	m_run_clock.Start();

    // Execute and signal Stepped
    m_simulator->TimeStep();

    m_run_clock.Stop();
}

void SimController::End()
{
	m_run_clock.Stop();
    PrintSummary();
    m_run_clock.Reset();
}
} // namespace stride
