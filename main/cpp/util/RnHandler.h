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
 *  Copyright 2019, Willem L, Kuylen E, Broeckhove J
 *  Copyright 2021, Willem L, Kuylen E, Libin P
 */

/**
 * @file
 * Header for the RnHandler class.
 */

#pragma once

#include <cmath>
#include <functional>

namespace stride {
namespace util {

/**
 * Draw random numbers [0,1] and perform binomial trials (or Bernouilli trials).
 */
class RnHandler
{
public:
        /// Constructor sets the transmission probability and random number generator.
        explicit RnHandler(std::function<double()> gen) : m_uniform01_generator(std::move(gen)) {}

        /// Make a draw on the uniform generator.
        double operator()() { return m_uniform01_generator(); }

        /// Perform binomial trial with given probability.
		bool Binomial(double probability_a)
		{
				return m_uniform01_generator() < probability_a;
		}

		/// Perform binomial trial with given the product of the given probabilities.
        bool Binomial(double probability_a, double probability_b)
        {
                return m_uniform01_generator() < probability_a * probability_b;
        }

private:
        /// Convert (exponential) rate into probability
        double RateToProbability(double rate) { return 1.0 - std::exp(-rate); }

private:
        std::function<double()> m_uniform01_generator; ///< Random number generator: double in [0.0, 1.0)
};

} // namespace util
} // namespace stride
