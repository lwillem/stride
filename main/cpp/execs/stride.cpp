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
 * Main program: command line handling.
 */

#include "SimController.h"
#include "util/FileSys.h"
#include "util/RunConfigManager.h"
#include "util/StringUtils.h"
#include "util/TimeStamp.h"

#include "util/Ptree.h"
#include <tclap/CmdLine.h>
#include <iostream>
#include <regex>
#include <string>
#include <vector>

using namespace std;
using namespace stride;
using namespace stride::util;
using namespace TCLAP;

/// Main program of the stride simulator.
int main(int argc, char** argv)
{
        int exitStatus = EXIT_SUCCESS;

        try {
                // -----------------------------------------------------------------------------------------
                // Parse command line (parameters displayed in --help in reverse order to order below).
                // -----------------------------------------------------------------------------------------
                CmdLine cmd("stride", ' ', "3.0");

                string sc = "Specifies the run configuration parameters. The format may be  is -c <file> ."
                            "\nDefaults to -c file=./config/run_default.xml";
                ValueArg<string> configArg("c", "config", sc, false, "config/run_default.xml", "CONFIGURATION", cmd);


                cmd.parse(argc, static_cast<const char* const*>(argv));

                // -----------------------------------------------------------------------------------------
                // Get configuration.
                // -----------------------------------------------------------------------------------------
                auto  config = configArg.getValue();
                ptree configPt = FileSys::ReadPtreeFile(config);;


                // -----------------------------------------------------------------------------------------
                // config and run simulation in cli
                // -----------------------------------------------------------------------------------------

                // add timestamp if no output prefix specified
				if (configPt.get<string>("run.output_prefix", "").empty()) {
						configPt.put("run.output_prefix", TimeStamp().ToTag().append("/"));
				}

                // sort the configuration details
				configPt.sort();

				// activate the controller
				SimController(configPt).Control();


        } catch (exception& e) {
                exitStatus = EXIT_FAILURE;
                cerr << "\nEXCEPTION THROWN: " << e.what() << endl;
        } catch (...) {
                exitStatus = EXIT_FAILURE;
                cerr << "\nEXCEPTION THROWN: Unknown exception." << endl;
        }
        return exitStatus;
}
