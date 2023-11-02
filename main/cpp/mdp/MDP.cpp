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
 * Implementation for the MDP class.
 */

#include "MDP.h"

#include "contact/ContactType.h"
#include "mdp/MDPRunner.h"
#include "pop/Person.h"
#include "pop/Population.h"
#include "sim/Sim.h"
#include "util/FileSys.h"
#include "util/RunConfigManager.h"
#include "util/TimeStamp.h"

#include <utility>

namespace stride {

using namespace std;
using namespace stride::util;
using namespace EventLogMode;
using namespace stride::ContactType;

MDP::MDP()
    : m_config(), m_simulator(nullptr), m_runner(nullptr), m_rnMan(), m_age_groups(), m_childless_age_groups(),
      m_mRNA_properties(nullptr), m_adeno_properties(nullptr)
{
}

MDP::~MDP() {
    cout << "Deleting MDP attributes..." << endl;
    if (!m_simulator->GetPopulation()->empty()) { ClearSimulation(); }
}

void MDP::Create(const std::string& configPath,
                 std::shared_ptr<VaccineProperties> mRNA_properties, std::shared_ptr<VaccineProperties> adeno_properties,
                 int seed, const std::string& outputDir, const std::string& outputPrefix, bool childless, double uptake)
{
    boost::property_tree::ptree configPt = FileSys::ReadPtreeFile(configPath);
    MDP::Create_(configPt, seed, outputDir, outputPrefix, childless, uptake);
    m_mRNA_properties = mRNA_properties;
    m_adeno_properties = adeno_properties;
}

void MDP::Create_(const boost::property_tree::ptree& config, int seed,
                  const std::string& outputDir, const std::string& outputPrefix, bool childless, double uptake)
{
    // Update the config
    m_config = config;
    // Update the output directory
    if (!outputDir.empty()) {
        m_config.put("run.output_prefix", outputDir);
    }
    // From SimController.Control
    // add timestamp if no output prefix specified
    else if (m_config.get<string>("run.output_prefix", "").empty()) {
        m_config.put("run.output_prefix", TimeStamp().ToTag().append("/"));
    }
    // Sort the configuration details
    m_config.sort();
    // Update the control helper since we didn't supply it at instantiation
    ControlHelper::m_config = m_config;
    ControlHelper::m_output_prefix = outputDir;

    // From SimController.Control
    // -----------------------------------------------------------------------------------------
    // Prelims.
    // -----------------------------------------------------------------------------------------
    CheckEnv();
    CheckOutputPrefix();
    InstallLogger();
    LogStartup();

    // Update the output directory to the prefix directory for current simulation
    if (!outputPrefix.empty()) {
        m_config.put("run.output_prefix", outputPrefix);
        // Update the control helper since we didn't supply it at instantiation
        ControlHelper::m_config = m_config;
        ControlHelper::m_output_prefix = outputPrefix;
    }

    if (seed != NULL) {
        m_stride_logger->info("Setting seed {}", seed);
        m_config.put("run.rng_seed", seed);
    }

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
    // Sim scenario: step , build a runner, register viewers.
    // -----------------------------------------------------------------------------------------
    auto runner = make_shared<MDPRunner>(m_config, m_simulator);
    RegisterViewers(runner);
    m_runner = runner;

    // -----------------------------------------------------------------------------------------
    // Vaccines: Create the age groups for vaccine sampling later
    // -----------------------------------------------------------------------------------------
    m_rnMan = rnMan;
    if (uptake != 1) { CreateHouseholdMapping(uptake); }
    else if (!childless) { CreateAgeGroups(); }
    else { CreateChildlessAgeGroups(); }
}

void MDP::UpdateCntReduction(std::vector<double> workplace_distancing, std::vector<double> community_distancing,
                             std::vector<double> collectivity_distancing)
{
    m_simulator->GetCalendar()->UpdateCntReduction(workplace_distancing, community_distancing, collectivity_distancing);
}

unsigned int MDP::Simulate(unsigned int numDays)
{
    // Run the simulation the given number of days
    for (unsigned int i = 0; i < numDays; i++) {
        m_runner->Step();
    }
    // Return the number of infected
    return m_simulator->GetPopulation()->GetTotalInfected();
}

unsigned int MDP::SimulateDay()
{
    // Run the simulation for a day
    return MDP::Simulate(1);
}

unsigned int MDP::SimulateVaccinate(unsigned int numDays, unsigned int availableVaccines,
                                    AgeGroup ageGroup, VaccineType vaccineType)
{
    // Run the simulation the given number of days
    for (unsigned int i = 0; i < numDays; i++) {
        // Vaccinate
        Vaccinate(availableVaccines, ageGroup, vaccineType);
        // Simulate
        m_runner->Step();
    }
    // Return the number of infected
    return m_simulator->GetPopulation()->GetTotalInfected();
}

void MDP::Vaccinate(unsigned int availableVaccines, AgeGroup ageGroup, VaccineType vaccineType)
{
    // Get people availableVaccines people within the given age group
    std::shared_ptr<Population> pop = m_simulator->GetPopulation();
    // Sample from the given age group (none of these are vaccinated yet)
    std::vector<unsigned int> sampled = SampleAgeGroup(ageGroup, availableVaccines);
    // Get the vaccine properties for the given vaccine type
    std::shared_ptr<VaccineProperties> vaccineProperties;
    // mRNA
    if (vaccineType == mRNA) { vaccineProperties = m_mRNA_properties; }
    // Adeno
    else if (vaccineType == adeno) { vaccineProperties = m_adeno_properties; }

    // Vaccinate the sampled people
    for (unsigned int id : sampled) {
        Person& p = pop->at(id);
        auto vaccine = vaccineProperties->GetVaccine();
        p.SetVaccine(vaccine);
        // Another person vaccinated from the age group
        m_vaccinated_age_groups[ageGroup] += 1;
    }
    // Log info
    m_stride_logger->info("[Vaccinate] {}/{} people vaccinated from age group {} with vaccine {}",
                          sampled.size(), availableVaccines, ageGroup, vaccineType);
}

void MDP::VaccinateChildless(unsigned int availableVaccines, ChildlessAgeGroup ageGroup, VaccineType vaccineType)
{
    // Get people availableVaccines people within the given age group
    std::shared_ptr<Population> pop = m_simulator->GetPopulation();
    // Sample from the given age group (none of these are vaccinated yet)
    std::vector<unsigned int> sampled = SampleAgeGroup(ageGroup, availableVaccines);
    // Get the vaccine properties for the given vaccine type
    std::shared_ptr<VaccineProperties> vaccineProperties;
    // mRNA
    if (vaccineType == mRNA) { vaccineProperties = m_mRNA_properties; }
        // Adeno
    else if (vaccineType == adeno) { vaccineProperties = m_adeno_properties; }

    // Vaccinate the sampled people
    for (unsigned int id : sampled) {
        Person& p = pop->at(id);
        auto vaccine = vaccineProperties->GetVaccine();
        p.SetVaccine(vaccine);
        // Another person vaccinated from the age group
        m_vaccinated_childless_age_groups[ageGroup] += 1;
    }
    // Log info
    m_stride_logger->info("[VaccinateChildless] {}/{} people vaccinated from age group {} with vaccine {}",
                          sampled.size(), availableVaccines, ageGroup, vaccineType);
}

void MDP::End()
{
    m_runner->End();
}

unsigned int MDP::GetNumberOfDays()
{
    const auto numDays = m_config.get<unsigned int>("run.num_days");
    return numDays;
}

unsigned int MDP::GetPopulationSize()
{
    return m_simulator->GetPopulation()->size();
}

std::map<AgeGroup, unsigned int> MDP::GetAgeGroupSizes()
{
    // Note: takes the number of people per age group that are not yet vaccinated
    //      => these values may change over time as people get vaccinated
    std::map<AgeGroup, unsigned int> groupCounts;
    for (AgeGroup ageGroup : AllAgeGroups) { groupCounts[ageGroup] = m_age_groups[ageGroup].size(); }
    return groupCounts;
}

std::map<ChildlessAgeGroup, unsigned int> MDP::GetChildlessAgeGroupSizes()
{
    // Note: takes the number of people per age group that are not yet vaccinated
    //      => these values may change over time as people get vaccinated
    std::map<ChildlessAgeGroup, unsigned int> groupCounts;
    for (ChildlessAgeGroup ageGroup : AllChildlessAgeGroups) { groupCounts[ageGroup] = m_childless_age_groups[ageGroup].size(); }
    return groupCounts;
}

unsigned int MDP::GetTotalInfected() const
{
    return m_simulator->GetPopulation()->GetTotalInfected();
}

unsigned int MDP::CountInfectedCases() const
{
    return m_simulator->GetPopulation()->CountInfectedCases();
}

unsigned int MDP::CountExposedCases() const
{
    return m_simulator->GetPopulation()->CountExposedCases();
}

unsigned int MDP::CountInfectiousCases() const
{
    return m_simulator->GetPopulation()->CountInfectiousCases();
}

unsigned int MDP::CountSymptomaticCases() const
{
    return m_simulator->GetPopulation()->CountSymptomaticCases();
}

unsigned int MDP::CountHospitalisedCases() const
{
    return m_simulator->GetPopulation()->CountHospitalisedCases();
}

unsigned int MDP::GetTotalHospitalised() const
{
    return m_simulator->GetPopulation()->GetTotalHospitalised();
}

unsigned int MDP::GetAtRisk() const
{
    return m_simulator->GetPopulation()->GetAtRisk();
}

void MDP::CreateAgeGroups()
{
    m_stride_logger->info("Creating age groups...");
    // Create an empty mapping for each age group
    for (AgeGroup ageGroup : AllAgeGroups) { m_age_groups[ageGroup] = std::vector<unsigned int>(); }
    for (AgeGroup ageGroup : AllAgeGroups) { m_vaccinated_age_groups[ageGroup] = 0; }
    // Iterate over the population and add people to their age group
    std::shared_ptr<Population> pop = m_simulator->GetPopulation();
    for (Person& p : *pop) { m_age_groups[GetAgeGroup(p.GetAge())].push_back(p.GetId()); }
    // Remove unused capacity from the vectors and shuffle the values
    for (AgeGroup ageGroup : AllAgeGroups) {
        m_age_groups[ageGroup].shrink_to_fit();
        m_rnMan.Shuffle(m_age_groups[ageGroup], 0U);
    }
}

void MDP::CreateChildlessAgeGroups()
{
    m_stride_logger->info("Creating childless age groups...");
    // Create an empty mapping for each age group
    for (ChildlessAgeGroup ageGroup : AllChildlessAgeGroups) { m_childless_age_groups[ageGroup] = std::vector<unsigned int>(); }
    for (ChildlessAgeGroup ageGroup : AllChildlessAgeGroups) { m_vaccinated_childless_age_groups[ageGroup] = 0; }
    // Iterate over the population and add people to their age group
    std::shared_ptr<Population> pop = m_simulator->GetPopulation();
    for (Person& p : *pop) {
        // Only store adults
        ChildlessAgeGroup ageGroup = GetChildlessAgeGroup(p.GetAge());
        if (ageGroup != ChildlessAgeGroup::children_c) { m_childless_age_groups[ageGroup].push_back(p.GetId()); }
    }
    // Remove unused capacity from the vectors and shuffle the values
    for (ChildlessAgeGroup ageGroup : AllChildlessAgeGroups) {
        m_childless_age_groups[ageGroup].shrink_to_fit();
        m_rnMan.Shuffle(m_childless_age_groups[ageGroup], 0U);
    }
}

std::vector<unsigned int> MDP::SampleAgeGroup(AgeGroup ageGroup, unsigned int samples)
{
    // Get the age group to sample from
    std::vector<unsigned int> sampled;
    // Take the last elements of the shuffled age group and remove them from the vector
    while (!m_age_groups[ageGroup].empty() && sampled.size() < samples) {
        sampled.push_back(m_age_groups[ageGroup].back());
        m_age_groups[ageGroup].pop_back();
    }
    return sampled;
}

std::vector<unsigned int> MDP::SampleAgeGroup(ChildlessAgeGroup ageGroup, unsigned int samples)
{
    // Get the age group to sample from
    std::vector<unsigned int> sampled;
    // Take the last elements of the shuffled age group and remove them from the vector
    while (!m_childless_age_groups[ageGroup].empty() && sampled.size() < samples) {
        sampled.push_back(m_childless_age_groups[ageGroup].back());
        m_childless_age_groups[ageGroup].pop_back();
    }
    return sampled;
}

void MDP::CreateHouseholdMapping(double uptake)
{
    m_stride_logger->info("Creating household mapping...");
    int fullPopSize = 0;
    std::unordered_map<unsigned int, unsigned int> householdSizes;
    std::vector<unsigned int> householdIds;
    // Iterate over contact pools
    auto& poolSys = m_simulator->GetPopulation()->RefPoolSys();
    const util::SegmentedVector<ContactPool>& householdPools = poolSys.CRefPools(ContactType::Id::Household);

    // Iterate over contact pools
    for (const ContactPool& h : householdPools) {
        const unsigned int hhId = h.GetId();
        const unsigned int hSize = h.size();
        householdSizes[hhId] = hSize;
        householdIds.push_back(hhId);
        fullPopSize += hSize;
    }
    // Remove unused capacity from container
    householdIds.shrink_to_fit();

    m_stride_logger->info("Sampling household update...");
    // Determine the approximate amount of persons from the population that need to be sampled
    unsigned int maxSampleSize = (int)std::round(uptake * fullPopSize);
    unsigned int numPools = householdIds.size();
    cout << "\tcontact pools: " << numPools << ", population size: " << fullPopSize << ", uptake: " << uptake <<
            " ==> requested sample size " << maxSampleSize << endl;
    std::vector<unsigned int> preSampledHouseholds (householdIds);
    m_rnMan.Shuffle(preSampledHouseholds, 0U);

    std::vector<unsigned int> samplePools;  // TODO remove once test cleared
    std::vector<unsigned int> sampleIds;
    unsigned int sampleSize = 0;
    unsigned int bestMatch = 0;
    unsigned int bestSize = 0;

    while (sampleSize < maxSampleSize) {
        // No more pools left to sample
        if (preSampledHouseholds.empty()) { break; }

        unsigned int newPoolId = preSampledHouseholds.back();
        preSampledHouseholds.pop_back();
        unsigned int poolSize = householdSizes[newPoolId];
        unsigned int newSize = sampleSize + poolSize;

        // Done
        if (newSize == maxSampleSize) {
            sampleSize = newSize;
            samplePools.push_back(poolSize);
            sampleIds.push_back(newPoolId);
        }
        // If size is larger than allowed, resample but keep best sample
        else if (newSize > maxSampleSize) {
            if (bestSize < newSize) { continue; }
            bestSize = newSize;
            bestMatch = newPoolId;
        }
        // Gather samples
        else {
            sampleSize = newSize;
            samplePools.push_back(poolSize);
            sampleIds.push_back(newPoolId);
        }
    }

    // Should normally not trigger for large populations with many contact pools
    if (bestSize != 0 and sampleSize != maxSampleSize) {
        unsigned int diffWithout = abs(int(maxSampleSize - sampleSize));
        unsigned int diffWith = abs(int(maxSampleSize - bestSize));
        cout << "\tCouldn't find a perfect match, deciding remaining best household... diffWithout = " <<
                diffWithout << ", diffWith = " << diffWith << endl;
        if (diffWith < diffWithout) {
            cout << "\tChoosing with remaining household" << endl;
            sampleSize = bestSize;
            samplePools.push_back(householdSizes[bestMatch]);
            sampleIds.push_back(bestMatch);
        }
    }
    cout << "\tSampled contact pools: " << sampleIds.size() << ", sampled persons: " << sampleSize << endl;

    // Create age groups, with only the people in the selected contact pools
    m_stride_logger->info("Creating age groups from households with uptake...");
    // Create an empty mapping for each age group
    for (AgeGroup ageGroup : AllAgeGroups) { m_age_groups[ageGroup] = std::vector<unsigned int>(); }
    for (AgeGroup ageGroup : AllAgeGroups) { m_vaccinated_age_groups[ageGroup] = 0; }
    // Iterate over the contact pools and add people to their age group
    for (const unsigned int hhId : sampleIds) {
        const ContactPool& h = householdPools.at(hhId);
        for (const Person *p : h.GetPool()) { m_age_groups[GetAgeGroup(p->GetAge())].push_back(p->GetId()); }
    }
    // Remove unused capacity from the vectors and shuffle the values
    for (AgeGroup ageGroup : AllAgeGroups) {
        m_age_groups[ageGroup].shrink_to_fit();
        m_rnMan.Shuffle(m_age_groups[ageGroup], 0U);
    }
}

void MDP::ClearSimulation()
{
    // Clear the ContactPoolSys
    cout << "\tClearing ContactPoolSys..." << endl;
    m_simulator->GetPopulation()->RefPoolSys().ClearContactPools();

    cout << "\tClearing population..." << endl;
    m_simulator->GetPopulation()->clear();

    cout << "\tClearing simulator and runner..." << endl;
    m_runner->End();
    m_runner.reset();
    m_simulator.reset();

    cout << "\tClearing age groups..." << endl;
    m_age_groups.clear();
    m_childless_age_groups.clear();

    cout << "\tClearing vaccine properties..." << endl;
    m_mRNA_properties.reset();
    m_adeno_properties.reset();
}

} // namespace stride
