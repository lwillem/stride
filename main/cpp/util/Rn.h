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

namespace stride {
namespace util {


class Rn {

  public:
	Rn() : engine_() {}

    trng::lcg64& engine()
     {
         return engine_;
     }

 	const trng::lcg64& engine() const
 	{
 			return engine_;
 	}

     template <typename Iter>
     void shuffle(Iter first, Iter last)
     {
         std::shuffle(first, last, engine_);
     }

  private:
       trng::lcg64 engine_;

 }; // end class

} // end namespace util
} // end namespace stride


