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
 * Implementation of the Gamma CDF and quantile functions.
 */

#include "GammaDistribution.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace stride {
namespace util {

namespace {

constexpr int    MAX_ITER = 200;
constexpr double EPS      = 1e-14;
const double     FPMIN    = std::numeric_limits<double>::min() / EPS;

/// Regularized lower incomplete gamma function P(a, x), via series expansion.
/// Accurate/efficient for x < a + 1.
double GammaSeries(double a, double x)
{
        if (x <= 0.0) {
                return 0.0;
        }
        double ap  = a;
        double sum = 1.0 / a;
        double del = sum;
        for (int n = 0; n < MAX_ITER; ++n) {
                ap += 1.0;
                del *= x / ap;
                sum += del;
                if (std::fabs(del) < std::fabs(sum) * EPS) {
                        break;
                }
        }
        return sum * std::exp(-x + a * std::log(x) - std::lgamma(a));
}

/// Regularized upper incomplete gamma function Q(a, x) = 1 - P(a, x), via a
/// continued fraction (modified Lentz's method). Accurate/efficient for x >= a + 1.
double GammaContinuedFraction(double a, double x)
{
        double b = x + 1.0 - a;
        double c = 1.0 / FPMIN;
        double d = 1.0 / b;
        double h = d;
        for (int i = 1; i < MAX_ITER; ++i) {
                const double an = -static_cast<double>(i) * (static_cast<double>(i) - a);
                b += 2.0;
                d = an * d + b;
                if (std::fabs(d) < FPMIN) {
                        d = FPMIN;
                }
                c = b + an / c;
                if (std::fabs(c) < FPMIN) {
                        c = FPMIN;
                }
                d                = 1.0 / d;
                const double del = d * c;
                h *= del;
                if (std::fabs(del - 1.0) < EPS) {
                        break;
                }
        }
        return std::exp(-x + a * std::log(x) - std::lgamma(a)) * h;
}

} // namespace

double GammaCdf(double x, double shape, double scale)
{
        if (x <= 0.0) {
                return 0.0;
        }
        const double y = x / scale;
        if (y < shape + 1.0) {
                return GammaSeries(shape, y);
        }
        return 1.0 - GammaContinuedFraction(shape, y);
}

double GammaQuantile(double p, double shape, double scale)
{
        if (p <= 0.0) {
                return 0.0;
        }
        if (p >= 1.0) {
                // unreachable in practice (the truncated-[0,1] sampling call site never
                // reaches p=1), but avoid std::numeric_limits<double>::infinity(), which
                // is UB under this project's -ffast-math build flags.
                return std::numeric_limits<double>::max();
        }

        // Bracket the root: the CDF is 0 at x=0 and increases monotonically to 1.
        double lo = 0.0;
        double hi = std::max(shape * scale, 1.0);
        while (GammaCdf(hi, shape, scale) < p) {
                hi *= 2.0;
        }

        // Safeguarded Newton-Raphson: falls back to a bisection step whenever the
        // Newton step would leave the current bracket, so this is guaranteed to
        // converge for a monotonic, smooth CDF like the gamma distribution's.
        double x = 0.5 * (lo + hi);
        for (int i = 0; i < MAX_ITER; ++i) {
                const double f = GammaCdf(x, shape, scale) - p;
                if (f < 0.0) {
                        lo = x;
                } else {
                        hi = x;
                }

                double xNext;
                if (x > 0.0) {
                        // PDF(x) = x^(shape-1) * exp(-x/scale) / (Gamma(shape) * scale^shape)
                        const double logPdf = (shape - 1.0) * std::log(x) - x / scale - std::lgamma(shape) - shape * std::log(scale);
                        const double pdf    = std::exp(logPdf);
                        xNext               = (pdf > 0.0) ? x - f / pdf : 0.5 * (lo + hi);
                } else {
                        xNext = 0.5 * (lo + hi);
                }
                if (!(xNext > lo && xNext < hi)) {
                        xNext = 0.5 * (lo + hi);
                }

                if (std::fabs(xNext - x) < EPS * std::fabs(xNext) || (hi - lo) < EPS) {
                        return xNext;
                }
                x = xNext;
        }
        return x;
}

} // namespace util
} // namespace stride
