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
 * Header for the InfectorMap.
 */

#pragma once

#include "contact/Infector.h"
#include "contact/InfectorExec.h"

#include <map>

#include "EventLogMode.h"

namespace stride {

class Calendar;
class Population;

/**
 * Mechanism to select the appropriate Infector template to execute.
 */
class InfectorMap : public std::map<stride::EventLogMode::Id, InfectorExec*>
{
public:
        /// Fully initialized.
        InfectorMap()
        {
                using namespace EventLogMode;

                this->emplace(Id::None, &Infector<Id::None>::Exec);
                this->emplace(Id::Incidence, &Infector<Id::Incidence>::Exec);
                this->emplace(Id::Transmissions, &Infector<Id::Transmissions>::Exec);
                this->emplace(Id::Participants, &Infector<Id::Participants>::Exec);
                this->emplace(Id::All, &Infector<Id::All>::Exec);
        }
};

} // namespace stride
