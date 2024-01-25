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
 * Interface of RnMan to manage random number generation (in parallel).
 */

#pragma once

#include <util/Rn.h>

//#include <trng/discrete_dist.hpp>
//#include <trng/lcg64.hpp>
//#include <trng/uniform01_dist.hpp>
//#include <trng/uniform_int_dist.hpp>
#include <functional>
#include <random>
#include <string>
#include <vector>

namespace stride {
namespace util {

/**
 * Manages random number generation in parallel (OpenMP) calculations.
 */
class RnMan : protected std::vector<util::Rn>
{
public:
        using ContainerType = std::vector<Rn>;
        using ContainerType::operator[];
        using ContainerType::at;
        using ContainerType::size;

public:
        /// Default constructor build empty manager.
        RnMan() : ContainerType(), m_seed_init(0U), m_stream_count(0U) {}

        /// Constructor.
		RnMan(unsigned long seed, const unsigned long stream_count)
			: ContainerType(stream_count),
			  m_seed_init(seed),
			  m_stream_count(stream_count)
		{
            Seed(m_seed_init);
		}

        /// No copying.
        RnMan(const RnMan&) = delete;

        /// No copy assignment.
        RnMan& operator=(const RnMan&) = delete;

        /// Return a generator for uniform doubles in [0, 1[ using i-th random engine.
        std::function<double()> GetUniform01Generator(unsigned int i = 0U)
        {
          	return std::bind(trng::uniform01_dist<double>(), std::ref(ContainerType::at(i).engine()));
        }

        /// Return a generator for uniform ints in [a, b[ (a < b) using i-th random engine.
        std::function<int()> GetUniformIntGenerator(int a, int b, unsigned int i = 0U)
        {
           	return std::bind(trng::uniform_int_dist(a, b), std::ref(ContainerType::at(i).engine()));
        }

        /// Return a generator for doubles from a Gamma distribution with a given shape and scale
        std::function<double()> GetGammaGenerator(double shape, double scale, unsigned int i = 0U)
		{
        	return std::bind(std::gamma_distribution<double>(shape, scale), std::ref(ContainerType::at(i).engine()));
		}

        /// Return generator for integers [0, n-1[ with non-negative weights p_j (i=0,..,n-1) using i-th random engine.
        template<typename It>
        std::function<int()> GetDiscreteGenerator(It begin, It end, unsigned int i = 0U)
        {
        	return std::bind(trng::discrete_dist(begin, end), std::ref(ContainerType::at(i).engine()));
        }

        /// Is this een empty (i.e. non-initialized RnMan)?
        bool IsEmpty() const { return ContainerType::empty() || (m_stream_count == 0U); }

        /// Random shuffle of vector of unsigned int indices using i-th engine.
        void Shuffle(std::vector<unsigned int>& indices, unsigned int i)
        {
                ContainerType::at(i).shuffle(indices.begin(), indices.end());
        }

        /// Make weighted coin flip: <fraction> of the flips need to come up true.
        bool MakeWeightedCoinFlip(double fraction, unsigned int i = 0U);

private:

        /// Actual first-time seeding. Procedure varies according to engine type, see specialisations.
        void Seed(unsigned long rng_seed);

private:
        unsigned long  m_seed_init;     ///< Seed initializer used with RN engine.
        unsigned int   m_stream_count;  ///< Number of threads/streams set up with the engine.
};


} // namespace util
} // namespace stride
