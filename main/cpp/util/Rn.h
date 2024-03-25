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
 * Interface of Rn to manage one random number generator.
 */

#pragma once

#include <trng/discrete_dist.hpp>
#include <trng/lcg64.hpp>
#include <trng/uniform01_dist.hpp>
#include <trng/uniform_int_dist.hpp>

#include <random>

namespace stride {
namespace util {


class Rn {

  public:
	Rn() : m_engine() {}

	double SampleUniform01(){
		return uniform01(m_engine);
	}

	/// Perform binomial trial with given probability.
	bool Binomial(double probability_a)
	{
		return SampleUniform01() < probability_a;
	}

	/// Perform binomial trial with the product of the given probabilities.
    bool Binomial(double probability_a, double probability_b)
    {
    	return SampleUniform01() < probability_a * probability_b;
    }

    trng::lcg64& GetEngine()
    {
       return m_engine;
    }

 	const trng::lcg64& GetEngine() const
 	{
 		return m_engine;
 	}

     template <typename Iter>
     void shuffle(Iter first, Iter last)
     {
         std::shuffle(first, last, m_engine);
     }

     /// Return a generator function for uniform integers in [a, b[ (a < b)
	 std::function<int()> GetUniformIntGenerator(int a, int b)
	 {
	 	return std::bind(trng::uniform_int_dist(a, b), std::ref(GetEngine()));
	 }

	 /// Return a generator function for doubles from a Gamma distribution with a given shape and scale
	 std::function<double()> GetGammaGenerator(double shape, double scale)
	 {
		return std::bind(std::gamma_distribution<double>(shape, scale), std::ref(GetEngine()));
	 }

     /// Random shuffle of vector of unsigned integers indices
     void Shuffle(std::vector<unsigned int>& indices)
     {
     	shuffle(indices.begin(), indices.end());
     }

  private:

       trng::lcg64 m_engine;                     /// random number engine
       trng::uniform01_dist<double> uniform01;   /// uniform distribution between 0 and 1

 }; // end class

} // end namespace util
} // end namespace stride


