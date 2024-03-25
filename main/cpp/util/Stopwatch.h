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
 * Definition of Stopwatch.
 */

#pragma once

#include "TimeToString.h"
#include <chrono>
#include <string>

namespace stride {
namespace util {

/**
 * Provides a stopwatch interface to time: it accumulates time between
 * start/stop pairs.
 */
class Stopwatch
{
public:
        using TClock = std::chrono::steady_clock;

        /// Constructor initializes stopwatch.
        explicit Stopwatch(std::string name = "stopwatch", bool running = false)
            : m_accumulated(TClock::duration::zero()), m_last_start(), m_name(std::move(name)), m_running(running)
        {
                if (m_running) {
                        m_last_start = std::chrono::steady_clock::now();
                }
        }

        /// Starts stopwatch if it was stopped.
        Stopwatch& Start()
        {
                if (!m_running) {
                        m_running    = true;
                        m_last_start = TClock::now();
                }
                return *this;
        }

        /// Stops the stopwatch if it was running.
        Stopwatch& Stop()
        {
                if (m_running) {
                        m_accumulated += (TClock::now() - m_last_start);
                        m_running = false;
                }
                return *this;
        }

        /// Resets stopwatch i.e. stopwatch is stopped and time accumulator is cleared.
        Stopwatch& Reset()
        {
                m_accumulated = std::chrono::steady_clock::duration::zero();
                m_running     = false;
                return *this;
        }

        /// Reports whether stopwatch has been started.
        bool IsRunning() const { return (m_running); }

        /// Return name of this stopwatch
        std::string GetName() const { return m_name; }

        /// Returns the accumulated value without altering the stopwatch state.
        typename std::chrono::steady_clock::duration Get() const
        {
                auto fTemp = m_accumulated;
                if (m_running) {
                        fTemp += (TClock::now() - m_last_start);
                }
                return fTemp;
        }

        /// Returns integer representation of readout
		unsigned int ToUnsignedInteger() const
		{

				using namespace std::chrono;
				milliseconds dur    = duration_cast<milliseconds>(Get());
				return static_cast<unsigned int>(dur.count());

		}

        /// Returns string representation of readout
        std::string ToString() const
        {
                using namespace std;
                using namespace std::chrono;

                string                          colon_string;
                typedef typename TClock::period TPeriod;
                if (ratio_less_equal<TPeriod, micro>::value) {
                        microseconds d = duration_cast<microseconds>(Get());
                        colon_string   = TimeToString::ToColonString(d);
                } else if (ratio_less_equal<TPeriod, milli>::value) {
                        milliseconds d = duration_cast<milliseconds>(Get());
                        colon_string   = TimeToString::ToColonString(d);
                } else {
                        seconds d    = duration_cast<seconds>(Get());
                        colon_string = TimeToString::ToColonString(d);
                }
                return colon_string;
        }

private:
        typename TClock::duration   m_accumulated;
        typename TClock::time_point m_last_start;
        std::string            m_name;
        bool                   m_running;
};

} // namespace util
} // namespace stride
