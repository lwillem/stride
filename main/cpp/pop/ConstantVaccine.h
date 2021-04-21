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
 *  Copyright 2020, Libin P, Willem L
 */

/**
 * @file
 * Header for the ConstantVaccine class.
 */

#pragma once

namespace stride {


class ConstantVaccine : public Vaccine 
{
public:
        struct Properties {
            std::string m_id;
            double m_ve_susceptible;
            double m_ve_infectiousness;
            double m_ve_severe;
        };

        ///
        explicit ConstantVaccine(std::shared_ptr<Properties> properties) : m_properties(properties)
        {
        }

        ~ConstantVaccine() {}

public:
        ///
        virtual double GetVeSusceptible() const {
            return m_properties->m_ve_susceptible;
        }

        ///
        virtual double GetVeInfectiousness() const { 
            return m_properties->m_ve_infectiousness;
        }

        ///
        virtual double GetVeSevere() const { 
            return m_properties->m_ve_severe;
        }

private:
        std::shared_ptr<Properties> m_properties;
};

} // namespace stride
