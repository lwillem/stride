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
 * Implementation of RnMan.
 */

#include "util/RnMan.h"
#include "StringUtils.h"

#include <cctype>
#include <sstream>
#include <stdexcept>
#include <iostream>

using namespace std;

namespace stride {
namespace util {

bool RnMan::MakeWeightedCoinFlip(double fraction, unsigned int i)
{
        array<double, 2> weights{ 1.0 - fraction, fraction};
        // -> 0, return is false -> not part of the fraction
        // -> 1, return is true -> part of the fraction
        auto dist = GetDiscreteGenerator(weights.begin(), weights.end(), i);
        return static_cast<bool>(dist());
}

void RnMan::Seed(unsigned long rng_seed)
{
        for (size_t i = 0; i < m_stream_count; ++i) {
                (*this)[i].engine().seed(rng_seed);
                (*this)[i].engine().split(m_stream_count, i);
        }
}


} // namespace util
} // namespace stride
