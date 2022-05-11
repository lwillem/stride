import matplotlib.pyplot as plt
import numpy as np

from scipy.stats import gamma

from run_util import get_mean_non_truncated_gamma
from plots import save_figure

def main():
    overdispersion_parameters = [10, 1, 0.6, 0.4, 0.2]

    # Cartoon figure for truncated Gamma distribution
    display_scenario_names = [r"$\alpha_{i}$ = " + str(a) for a in overdispersion_parameters]

    mean_transmission_probability = 0.08

    plt.axvline(mean_transmission_probability, color="lightgrey")
    for overdispersion in overdispersion_parameters:
        mean_transmission_probability_corrected = get_mean_non_truncated_gamma(mean_transmission_probability, overdispersion)
        shape = overdispersion
        scale = mean_transmission_probability_corrected / shape
        x = np.linspace(0, 1, 500)
        cdf1 = gamma.cdf(1, a=shape, scale=scale)
        cdf0 = gamma.cdf(0, a=shape, scale=scale)
        plt.plot(x, gamma.pdf(x, a=shape, scale=scale) / (cdf1 - cdf0))
    plt.xlabel("Individual transmission probability")
    plt.ylabel("Probability density")
    plt.legend(["mean"] + display_scenario_names)
    save_figure(".", "pdf_gamma_infectiousness")

    # Cartoon figure for Gamma distribution
    display_scenario_names = [r"$\alpha_{c}$ = " + str(a) for a in overdispersion_parameters]

    plt.axvline(1, color="lightgrey")
    for overdispersion in overdispersion_parameters:
        shape = overdispersion
        scale = 1 / shape
        x = np.linspace(0, 5, 500)
        plt.plot(x, gamma.pdf(x, a=shape, scale=scale))
    plt.xlabel("Individual contact factor")
    plt.ylabel("Probability density")
    plt.legend(["mean"] + display_scenario_names)

    save_figure(".", "pdf_gamma_contacts")

if __name__=="__main__":
    main()
