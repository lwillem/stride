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
 *  Copyright 2026
 */

/**
 * @file
 * CDF and quantile (inverse CDF) for the Gamma(shape, scale) distribution,
 * replacing boost::math::gamma_distribution for stride's one remaining use:
 * inverse-transform sampling from a gamma distribution truncated to [0, 1].
 * Implemented with the regularized incomplete gamma function (series
 * expansion + continued fraction) and a bracketed, safeguarded
 * Newton-Raphson root-finder, using only <cmath>.
 */

#pragma once

namespace stride {
namespace util {

/// Regularized incomplete gamma CDF: P(X <= x) for X ~ Gamma(shape, scale).
/// Requires shape > 0, scale > 0.
double GammaCdf(double x, double shape, double scale);

/// Inverse CDF (quantile) for Gamma(shape, scale): returns x such that
/// GammaCdf(x, shape, scale) == p. Requires shape > 0, scale > 0, 0 <= p <= 1.
double GammaQuantile(double p, double shape, double scale);

} // namespace util
} // namespace stride
