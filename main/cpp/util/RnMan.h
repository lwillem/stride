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
        using std::vector<Rn>::operator[];
        using std::vector<Rn>::at;
        using std::vector<Rn>::size;

public:
        /// Default constructor build empty manager.
        RnMan() : std::vector<Rn>() {}

        /// Constructor.
		RnMan(unsigned long rng_seed, const unsigned long stream_count)
			: std::vector<Rn>(stream_count)
		{
			// seed random number generator(s)
			for (size_t i = 0; i < size(); ++i) {
				(*this)[i].GetEngine().seed(rng_seed);
				(*this)[i].GetEngine().split(size(), i);
			}
		}

        /// No copying.
        RnMan(const RnMan&) = delete;

        /// No copy assignment.
        RnMan& operator=(const RnMan&) = delete;

};


} // namespace util
} // namespace stride
