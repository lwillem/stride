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
 *  Copyright 2026, Willem L
 */

/**
 * @file
 * A boost::property_tree::ptree-like wrapper around pugixml, covering the
 * subset of ptree's API actually used by stride: dotted-path get/put with
 * defaults, get_optional, get_child, iteration over direct children, sort,
 * and XML (de)serialization from/to string or file. Kept dependency-free
 * from Boost so config parsing no longer needs boost::property_tree.
 */

#pragma once

#include <pugixml.hpp>

#include <algorithm>
#include <cctype>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace stride {
namespace util {

class ptree
{
public:
        using value_type = std::pair<std::string, ptree>;

        ptree() = default;
        ptree(const ptree& other) { m_doc.reset(other.m_doc); }
        ptree(ptree&& other) noexcept = default;
        ptree& operator=(const ptree& other)
        {
                if (this != &other) {
                        m_doc.reset(other.m_doc);
                }
                return *this;
        }
        ptree& operator=(ptree&& other) noexcept = default;
        ~ptree()                                 = default;

        /// Get the value at path, converted to T. Throws std::runtime_error if path is not found.
        template <typename T>
        T get(const std::string& path) const
        {
                const pugi::xml_node node = FindNode(path);
                if (!node) {
                        throw std::runtime_error("stride::util::ptree::get> path not found: " + path);
                }
                return Convert<T>(node.text().as_string());
        }

        /// Get the value at path, converted to T, or defaultValue if path is not found
        /// or cannot be converted to T.
        template <typename T>
        T get(const std::string& path, const T& defaultValue) const
        {
                const pugi::xml_node node = FindNode(path);
                if (!node) {
                        return defaultValue;
                }
                try {
                        return Convert<T>(node.text().as_string());
                } catch (const std::exception&) {
                        return defaultValue;
                }
        }

        /// Get the value at path, converted to T, or an empty optional if path is not found
        /// or cannot be converted to T.
        template <typename T>
        std::optional<T> get_optional(const std::string& path) const
        {
                const pugi::xml_node node = FindNode(path);
                if (!node) {
                        return std::nullopt;
                }
                try {
                        return Convert<T>(node.text().as_string());
                } catch (const std::exception&) {
                        return std::nullopt;
                }
        }

        /// Get the subtree at path (as an independent copy). Throws std::runtime_error
        /// if path is not found.
        ptree get_child(const std::string& path) const;

        /// Set (creating any missing intermediate elements) the value at path.
        /// Overwrites any existing value at that path.
        template <typename T>
        void put(const std::string& path, const T& value)
        {
                PutString(path, ToString(value));
        }

        /// Sort the direct children alphabetically by tag name (not recursive),
        /// matching boost::property_tree::ptree::sort().
        void sort();

        /// Iterator over (tag name, subtree) pairs of the direct children.
        class const_iterator
        {
        public:
                const_iterator() = default;
                explicit const_iterator(pugi::xml_node node) : m_node(node) { SkipToElement(); }

                bool operator!=(const const_iterator& other) const { return m_node != other.m_node; }
                bool operator==(const const_iterator& other) const { return m_node == other.m_node; }

                const_iterator& operator++()
                {
                        m_node = m_node.next_sibling();
                        SkipToElement();
                        return *this;
                }

                value_type operator*() const { return {m_node.name(), ptree::FromNode(m_node)}; }

        private:
                void SkipToElement()
                {
                        while (m_node && m_node.type() != pugi::node_element) {
                                m_node = m_node.next_sibling();
                        }
                }

                pugi::xml_node m_node;
        };

        const_iterator begin() const { return const_iterator(m_doc.first_child()); }
        const_iterator end() const { return const_iterator(pugi::xml_node()); }

        /// Parse XML from a string. Throws std::runtime_error on a parse error.
        static ptree FromXmlString(const std::string& xml);

        /// Serialize to an XML string (pretty-printed, 8-space indent).
        std::string ToXmlString() const;

        /// Parse XML from a file. Throws std::runtime_error on a parse error.
        static ptree FromXmlFile(const std::string& path);

        /// Serialize to an XML file (pretty-printed, 8-space indent).
        /// Throws std::runtime_error if the file cannot be written.
        void ToXmlFile(const std::string& path) const;

private:
        static ptree FromNode(const pugi::xml_node& node);
        pugi::xml_node FindNode(const std::string& path) const;
        void PutString(const std::string& path, const std::string& value);
        static bool ConvertBool(const std::string& s);

        template <typename T>
        static T Convert(const std::string& s)
        {
                if constexpr (std::is_same_v<T, std::string>) {
                        return s;
                } else if constexpr (std::is_same_v<T, bool>) {
                        return ConvertBool(s);
                } else {
                        std::istringstream iss(s);
                        T                  value{};
                        iss >> value;
                        if (iss.fail()) {
                                throw std::runtime_error("stride::util::ptree> cannot convert '" + s +
                                                          "' to requested type");
                        }
                        return value;
                }
        }

        template <typename T>
        static std::string ToString(const T& value)
        {
                if constexpr (std::is_same_v<T, std::string>) {
                        return value;
                } else if constexpr (std::is_same_v<T, bool>) {
                        return value ? "true" : "false";
                } else if constexpr (std::is_convertible_v<T, const char*>) {
                        return std::string(value);
                } else {
                        std::ostringstream oss;
                        oss << value;
                        return oss.str();
                }
        }

        pugi::xml_document m_doc;
};

} // namespace util
} // namespace stride
