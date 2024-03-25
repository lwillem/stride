# noinspection PyUnresolvedReferences
import pylibstride as stride
import random


if __name__ == '__main__':

    # The simulation file to load
    # TODO: don't forget to change the paths to the data files to use!
    default_file = "py_test/run_default.xml"

    # Create an instance of the MDP
    mdp = stride.MDP()
    mdp.Create(default_file)  # TODO: use same seed in python and c++?

    # Get the number of days to run the experiment for (from the configuration file)
    numDays = mdp.GetNumberOfDays()
    # The step size in days between the next action the bandit can take
    step_size = 10

    print(f"All age groups {stride.AllAgeGroups}")

    def random_vaccinate():
        group = random.choice(stride.AllAgeGroups)
        v_type = random.choice(stride.AllVaccineTypes[1:])  # ignore noVaccine
        return group, v_type

    # Run the simulation and vaccinate people
    for i in range(numDays // step_size):
        # Bandit should select action to vaccinate instead of random choice
        age_group, vaccine_type = random_vaccinate()
        print(f"Vaccinating {age_group.name} with {vaccine_type.name}")

        # Execute chosen vaccination for step_size days
        for _ in range(step_size):
            # Vaccinate
            mdp.Vaccinate(availableVaccines=10000, ageGroup=age_group, vaccineType=vaccine_type)
            # Run the simulation
            infected = mdp.SimulateDay()

        print("Number infected:", infected)

    # Stop the mdp (simulator)
    mdp.End()
