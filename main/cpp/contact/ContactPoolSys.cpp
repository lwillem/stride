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
 *  Copyright 2017, 2018, Kuylen E, Willem L, Broeckhove J
 */

/**
 * @file
 * Core Population class.
 */

#include <boost/property_tree/ptree.hpp>

#include "ContactPoolSys.h"

using namespace std;
using namespace stride::ContactType;

namespace stride {

ContactPoolSys::ContactPoolSys(const boost::property_tree::ptree& ventilation) : m_currentContactPoolId(), m_sys(), m_ventilation()
{

auto ventilationHousehold = ventilation.get<double>("ventilation_reduction.Household",0);
auto ventilationK12School = ventilation.get<double>("ventilation_reduction.K12School",0);
auto ventilationCollege = ventilation.get<double>("ventilation_reduction.College",0);
auto ventilationWorkplace = ventilation.get<double>("ventilation_reduction.Workplace",0);
auto ventilationPrimaryCommunity = ventilation.get<double>("ventilation_reduction.PrimaryCommunity",0);
auto ventilationSecondaryCommunity = ventilation.get<double>("ventilation_reduction.SecondaryCommunity",0);
auto ventilationHouseholdCluster = ventilation.get<double>("ventilation_reduction.HouseholdCluster",0);
auto ventilationCollectivity = ventilation.get<double>("ventilation_reduction.Collectivity",0);

        for (Id typ : IdList) {
                m_sys[typ].emplace_back(ContactPool(0U, typ, 0U));
                m_currentContactPoolId[typ] = 1;
                if (typ == ContactType::Id::Household) {m_ventilation[typ] = ventilationHousehold;}
                else if (typ == ContactType::Id::K12School) {m_ventilation[typ] = ventilationK12School;}
                else if (typ == ContactType::Id::College) {m_ventilation[typ] = ventilationCollege;}
                else if (typ == ContactType::Id::Workplace) {m_ventilation[typ] = ventilationWorkplace;}
                else if (typ == ContactType::Id::PrimaryCommunity) {m_ventilation[typ] = ventilationPrimaryCommunity;}
                else if (typ == ContactType::Id::SecondaryCommunity) {m_ventilation[typ] = ventilationSecondaryCommunity;}
                else if (typ == ContactType::Id::HouseholdCluster) {m_ventilation[typ] = ventilationHouseholdCluster;}
                else if (typ == ContactType::Id::Collectivity) {m_ventilation[typ] = ventilationCollectivity;}

        }

}

ContactPool* ContactPoolSys::CreateContactPool(ContactType::Id typeId)
{
        return m_sys[typeId].emplace_back(m_currentContactPoolId[typeId]++, typeId, m_ventilation[typeId]);
}

} // namespace stride
