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
 *  Copyright 2017, Kuylen E, Willem L, Broeckhove J
 */

/**
 * @file.
 * Miscellaneous string utilities.
 */

#pragma once

#include <algorithm>
#include <cctype>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <regex>

#include <vector>

namespace stride {
namespace util {

/// All characters in string are digits (or string is empty)
inline bool CheckAllDigits(const std::string& s)
{
        bool status = true;
        for (const auto& e : s) {
                status = (status && isdigit(e));
                if (!status)
                        break;
        }
        return status;
}

/// Builds a value of type T representation from a string.
template <typename T>
inline T FromString(const std::string& s)
{
        std::stringstream ss(s);
        T                 t;
        ss >> t;
        return t;
}

/// Split a string (in order of occurence) by splitting it on the given delimiters.
/// Consecutive delimiters produce empty tokens (not compressed/skipped).
inline std::vector<std::string> Split(const std::string& s, const std::string& delimiters)
{
        std::vector<std::string> tokens;
        std::string::size_type   start = 0;
        std::string::size_type   pos;
        while ((pos = s.find_first_of(delimiters, start)) != std::string::npos) {
                tokens.push_back(s.substr(start, pos - start));
                start = pos + 1;
        }
        tokens.push_back(s.substr(start));
        return tokens;
}

/// Convert a string to T, throwing std::invalid_argument if the whole string
/// can't be consumed as a valid T.
template <typename T>
inline T LexicalCast(const std::string& s)
{
        std::istringstream iss(s);
        T                  value{};
        iss >> value;
        if (iss.fail() || !iss.eof()) {
                throw std::invalid_argument("stride::util::LexicalCast> cannot convert '" + s + "' to requested type");
        }
        return value;
}

template <>
inline std::string LexicalCast<std::string>(const std::string& s)
{
        return s;
}

/// Tokenize a string (in order of occurence) with the given delimiters.
/// Multiple consecutive delimiters do NOT define "empty" tokens; they are skipped.
/// @throws std::invalid_argument if the token val can't be converted to T
template <typename T>
inline std::vector<T> Tokenize(const std::string& str, const std::string& delimiters)
{
        std::vector<T> tokens;

        // Skip delimiters at beginning.
        std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
        // Find first non-delimiter.
        std::string::size_type pos = str.find_first_of(delimiters, lastPos);

        while (std::string::npos != pos || std::string::npos != lastPos) {
                // Found a token, add it to the vector.
                tokens.push_back(LexicalCast<T>(str.substr(lastPos, pos - lastPos)));
                // Skip delimiters.
                lastPos = str.find_first_not_of(delimiters, pos);
                // Find next non-delimiter.
                pos = str.find_first_of(delimiters, lastPos);
        }

        return tokens;
}

/// Builds a string representation of a value of type T.
template <typename T>
inline std::string ToString(const T& value)
{
        std::stringstream ss;
        ss << value;
        return ss.str();
}

template <>
inline std::string ToString<std::string>(const std::string& value)
{
        return value;
}

/// Stringify values (that are not strings) in a range and put them in a vector.
template <typename It>
inline std::vector<std::string> ToString(
    typename std::enable_if<!std::is_same<typename It::value_type, std::string>::value, It>::type first, It last)
{
        std::vector<std::string> v;
        for (It it = first; it < last; ++it) {
                v.emplace_back(ToString(*it));
        }
        return v;
}

/// Stringify values (that are strings - so no-op) in a range and put them in a vector.
template <typename It>
inline std::vector<std::string> ToString(
    typename std::enable_if<std::is_same<typename It::value_type, std::string>::value, It>::type first, It last)
{
        std::vector<std::string> v;
        std::copy(first, last, back_inserter(v));
        return v;
}

/// Builds a string representation with minimum width of a value of type T.
template <typename T>
inline std::string ToString(const T& value, int width, char fill = ' ')
{
        std::stringstream ss;
        ss << std::setw(width) << std::setfill(fill) << value;
        return ss.str();
}

/// Builds a string with lower case characters only.
inline std::string ToLower(const std::string& source)
{
        auto        lower = [](int c) -> int { return std::tolower(c); };
        std::string copy;
        std::transform(source.begin(), source.end(), std::back_inserter(copy), lower);
        return copy;
}

/// Builds a string with upper case characters only.
inline std::string ToUpper(const std::string& source)
{
        auto        upper = [](int c) -> int { return std::toupper(c); };
        std::string copy;
        std::transform(source.begin(), source.end(), std::back_inserter(copy), upper);
        return copy;
}

/// Trim characters at right end of string.
inline std::string TrimRight(const std::string& source, const std::string& t = " ")
{
        std::string str = source;
        return str.erase(str.find_last_not_of(t) + 1);
}

/// Trim characters at left end of string.
inline std::string TrimLeft(const std::string& source, const std::string& t = " ")
{
        std::string str = source;
        return str.erase(0, source.find_first_not_of(t));
}

/// Trim characters at both ends of string.
inline std::string Trim(const std::string& source, const std::string& t = " ")
{
        return TrimLeft(TrimRight(source, t), t);
}

template <typename T>
inline std::string intToDottedString(const T& value)
{
        std::string valueStr = std::to_string(value);

        std::string res;
        std::size_t rest = valueStr.length() % 3;

        res += valueStr.substr(0, rest);

        for (size_t i = rest; i < valueStr.length(); i += 3) {
                res += "." + valueStr.substr(i, 3);
        }

        if (res[0] == '.') {
                return res.substr(1);
        }

        return res;
}

inline void Replace(std::string& s, const std::string& pattern, const std::string& replace)
{
        if (pattern.empty()) {
                return;
        }
        std::string::size_type pos = 0;
        while ((pos = s.find(pattern, pos)) != std::string::npos) {
                s.replace(pos, pattern.length(), replace);
                pos += replace.length();
        }
}

inline bool IsSubstring(std::string& s, std::string& pattern){

	std::regex e(pattern);
	return(std::regex_search(s,e));

}

} // namespace util
} // namespace stride
