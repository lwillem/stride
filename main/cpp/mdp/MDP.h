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
 * Header for the MDP class.
 */

#pragma once

#include "contact/AgeContactProfiles.h"
#include "contact/EventLogMode.h"
#include "contact/InfectorExec.h"
#include "disease/PublicHealthAgency.h"
#include "disease/TransmissionProfile.h"
#include "disease/UniversalTesting.h"
#include "mdp/AgeGroup.h"
#include "mdp/Vaccines.h"
#include "util/RnMan.h"
#include "util/RnHandler.h"
#include "execs/ControlHelper.h"

#include <boost/property_tree/ptree.hpp>
#include <string>
#include <map>
#include <unordered_map>
#include <tuple>

namespace stride {

class Calendar;
class Population;
class Sim;
class MDPRunner;

/**
 * Markov Decision Process for the stride simulation.
 * Exposes the simulation to distribute vaccines to age groups of the population.
 */
class MDP : protected ControlHelper
{
public:
        /// Constructor for empty MDP.
        explicit MDP();

        ~MDP();

        /// Create an MDP (and the underlying simulation) from a given configuration
        void Create(const std::string& configPath,
                    std::shared_ptr<VaccineProperties> mRNA_properties,
                    std::shared_ptr<VaccineProperties> adeno_properties,
                    int seed = 0, const std::string& outputDir = "", const std::string& outputPrefix = "",
                    bool childless = false, double uptake = 1);

        /// Update the contact reduction vectors
        void UpdateCntReduction(std::vector<double> workplace_distancing, std::vector<double> community_distancing,
                                std::vector<double> collectivity_distancing);

        /// Simulate a given number of days in the simulation
        unsigned int Simulate(unsigned int numDays);

        /// Simulate a single day in the simulation
        unsigned int SimulateDay();

        /// Simulate multiple days of vaccinations
        unsigned int SimulateVaccinate(unsigned int numDays, unsigned int availableVaccines,
                                       AgeGroup ageGroup, VaccineType vaccineType);

        /// Vaccinate a given age group with the given vaccine type and number of available vaccines
        void Vaccinate(unsigned int availableVaccines, AgeGroup ageGroup, VaccineType vaccineType);
        void VaccinateChildless(unsigned int availableVaccines, ChildlessAgeGroup ageGroup, VaccineType vaccineType);

        /// Notify the simulator should stop
        void End();

        /// Get the number of days specified to runt the simulation for
        unsigned int GetNumberOfDays();

        /// Get the population size
        unsigned int GetPopulationSize();

        /// Get the size of the age groups
        std::map<AgeGroup, unsigned int> GetAgeGroupSizes();
        std::map<ChildlessAgeGroup, unsigned int> GetChildlessAgeGroupSizes();

        /// Get the number of vaccinated individuals per age group
        std::map<AgeGroup, unsigned int> GetVaccinatedAgeGroups() { return m_vaccinated_age_groups; }
        std::map<ChildlessAgeGroup, unsigned int> GetVaccinatedChildlessAgeGroups() { return m_vaccinated_childless_age_groups; }

        /// Get the cumulative number of cases.
        unsigned int GetTotalInfected() const;

        /// Get the current number of infected cases.
        unsigned int CountInfectedCases() const;

        /// Get the current number of exposed cases.
        unsigned int CountExposedCases() const;

        /// Get the current number of infectious cases.
        unsigned int CountInfectiousCases() const;

        /// Get the current number of symptomatic cases.
        unsigned int CountSymptomaticCases() const;

        /// Get the current number of hospitalised cases.
        unsigned int CountHospitalisedCases() const;

        /// Get the cumulative number of hospitalisations.
        unsigned int GetTotalHospitalised() const;

        /// Get the number of people at risk in the population
        unsigned int GetAtRisk() const;

        /// Memory management
        void ClearSimulation();

        /// TODO: get MDP state
        // MDPState GetState();

private:
        /// Create an MDP (and the underlying simulation) from a given configuration
        void Create_(const boost::property_tree::ptree& config, int seed,
                     const std::string& outputDir, const std::string& outputPrefix, bool childless, double uptake);
        /// Create a mapping of the age groups with the person IDs of people corresponding to those age groups.
        void CreateAgeGroups();
        void CreateChildlessAgeGroups();
        /// Sample IDs for a given age group
        std::vector<unsigned int> SampleAgeGroup(AgeGroup ageGroup, unsigned int samples);
        std::vector<unsigned int> SampleAgeGroup(ChildlessAgeGroup ageGroup, unsigned int samples);
        /// Create a mapping of all households to their persons, to initialise age group mappings
        void CreateHouseholdMapping(double uptake);

private:
        boost::property_tree::ptree m_config;                       ///< Configuration property tree
        std::shared_ptr<Sim> m_simulator;                           ///< The simulation
        std::shared_ptr<MDPRunner> m_runner;                        ///< The runner for the simulation
        util::RnMan m_rnMan;                                        ///< The random number manager
        std::map<AgeGroup, std::vector<unsigned int>> m_age_groups; ///< The IDs of people belonging to different age groups.
        std::map<ChildlessAgeGroup, std::vector<unsigned int>> m_childless_age_groups;
        std::map<AgeGroup, unsigned int> m_vaccinated_age_groups;   ///< The number of vaccinated individuals per age group.
        std::map<ChildlessAgeGroup, unsigned int> m_vaccinated_childless_age_groups;
        std::shared_ptr<VaccineProperties> m_mRNA_properties;       ///< The vaccine properties of the mRNA vaccine
        std::shared_ptr<VaccineProperties> m_adeno_properties;      ///< The vaccine properties of the adeno vaccine
};

} // namespace stride
