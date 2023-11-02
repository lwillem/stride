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
 * Header for the mRNA and Adeno Vaccines.
 */

#pragma once

#include "pop/Vaccine.h"
#include "pop/ConstantVaccine.h"
#include "pop/LinearVaccine.h"

#include <string>
#include <vector>
#include <memory>


namespace stride {

// Python does not need access to vaccine class and its methods,
//  only specifies which vaccine to generate for simplicity
/// Vaccines enumeration
enum VaccineType { noVaccine, mRNA, adeno };

/// All the vaccine types
static const std::vector<VaccineType> AllVaccineTypes = { noVaccine, mRNA, adeno };


/// Vaccine properties
struct VaccineProperties {
    virtual ~VaccineProperties() = default;

    virtual std::unique_ptr<Vaccine> GetVaccine() = 0;
};

// Constant
struct ConstantVaccineProperties : VaccineProperties {
    ConstantVaccineProperties(std::string& v_id, double ve_susceptible, double ve_infectiousness, double ve_severe) :
            id(v_id), ve_susceptible(ve_susceptible), ve_infectiousness(ve_infectiousness), ve_severe(ve_severe) { }

    std::string id;
    double ve_susceptible;
    double ve_infectiousness;
    double ve_severe;

    virtual std::unique_ptr<Vaccine> GetVaccine() {
        // std::printf("Creating Constant Vaccine\n");
        std::shared_ptr<ConstantVaccine::Properties> properties;
        properties = std::shared_ptr<ConstantVaccine::Properties>(
                        new ConstantVaccine::Properties{id, ve_susceptible, ve_infectiousness, ve_severe});
        // Create and return the vaccine
        return std::unique_ptr<ConstantVaccine>(new ConstantVaccine(properties));
    };
};

// Linear
struct LinearVaccineProperties : VaccineProperties {
    LinearVaccineProperties(std::string& v_id, double ve_susceptible, double ve_infectiousness, double ve_severe,
                            unsigned short max_ve_day) :
                            id(v_id), ve_susceptible(ve_susceptible), ve_infectiousness(ve_infectiousness),
                            ve_severe(ve_severe), max_ve_day(max_ve_day) { }

    std::string id;
    double ve_susceptible;
    double ve_infectiousness;
    double ve_severe;
    unsigned short max_ve_day;

    virtual std::unique_ptr<Vaccine> GetVaccine() {
        // std::printf("Creating Linear Vaccine\n");
        std::shared_ptr<LinearVaccine::Properties> properties;
        properties = std::shared_ptr<LinearVaccine::Properties>(
                new LinearVaccine::Properties{id, ve_susceptible, ve_infectiousness, ve_severe, max_ve_day});
        // Create and return the vaccine
        return std::unique_ptr<LinearVaccine>(new LinearVaccine(properties));
    };
};

} // namespace stride
