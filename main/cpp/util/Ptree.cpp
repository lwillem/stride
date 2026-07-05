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
 * Implementation of the ptree wrapper around pugixml.
 */

#include "Ptree.h"

namespace stride {
namespace util {

namespace {

constexpr unsigned int g_parse_flags  = pugi::parse_default | pugi::parse_trim_pcdata;
const char*            g_indent       = "        "; // 8 spaces, matches the previous boost::property_tree settings

std::vector<std::string> SplitPath(const std::string& path)
{
        std::vector<std::string> parts;
        if (path.empty()) {
                return parts;
        }
        std::string::size_type start = 0;
        while (true) {
                const auto dot = path.find('.', start);
                if (dot == std::string::npos) {
                        parts.push_back(path.substr(start));
                        break;
                }
                parts.push_back(path.substr(start, dot - start));
                start = dot + 1;
        }
        return parts;
}

} // namespace

pugi::xml_node ptree::FindNode(const std::string& path) const
{
        pugi::xml_node node = m_doc;
        for (const auto& part : SplitPath(path)) {
                node = node.child(part.c_str());
                if (!node) {
                        return pugi::xml_node();
                }
        }
        return node;
}

ptree ptree::FromNode(const pugi::xml_node& node)
{
        ptree result;
        for (pugi::xml_node child = node.first_child(); child; child = child.next_sibling()) {
                result.m_doc.append_copy(child);
        }
        return result;
}

ptree ptree::get_child(const std::string& path) const
{
        const pugi::xml_node node = FindNode(path);
        if (!node) {
                throw std::runtime_error("stride::util::ptree::get_child> path not found: " + path);
        }
        return FromNode(node);
}

void ptree::PutString(const std::string& path, const std::string& value)
{
        pugi::xml_node node = m_doc;
        for (const auto& part : SplitPath(path)) {
                pugi::xml_node child = node.child(part.c_str());
                if (!child) {
                        child = node.append_child(part.c_str());
                }
                node = child;
        }
        node.text().set(value.c_str());
}

bool ptree::ConvertBool(const std::string& s)
{
        std::string lower = s;
        std::transform(lower.begin(), lower.end(), lower.begin(),
                        [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        if (lower == "true" || lower == "1") {
                return true;
        }
        if (lower == "false" || lower == "0") {
                return false;
        }
        throw std::runtime_error("stride::util::ptree> cannot convert '" + s + "' to bool");
}

void ptree::sort()
{
        std::vector<pugi::xml_node> children;
        for (pugi::xml_node child = m_doc.first_child(); child; child = child.next_sibling()) {
                children.push_back(child);
        }
        std::sort(children.begin(), children.end(), [](const pugi::xml_node& a, const pugi::xml_node& b) {
                return std::string(a.name()) < std::string(b.name());
        });

        pugi::xml_document sorted;
        for (const auto& child : children) {
                sorted.append_copy(child);
        }
        m_doc.reset(sorted);
}

ptree ptree::FromXmlString(const std::string& xml)
{
        ptree                        result;
        const pugi::xml_parse_result parseResult = result.m_doc.load_string(xml.c_str(), g_parse_flags);
        if (!parseResult) {
                throw std::runtime_error(std::string("stride::util::ptree::FromXmlString> ") +
                                          parseResult.description());
        }
        return result;
}

std::string ptree::ToXmlString() const
{
        std::ostringstream oss;
        m_doc.save(oss, g_indent, pugi::format_default);
        return oss.str();
}

ptree ptree::FromXmlFile(const std::string& path)
{
        ptree                        result;
        const pugi::xml_parse_result parseResult = result.m_doc.load_file(path.c_str(), g_parse_flags);
        if (!parseResult) {
                throw std::runtime_error(std::string("stride::util::ptree::FromXmlFile> ") +
                                          parseResult.description());
        }
        return result;
}

void ptree::ToXmlFile(const std::string& path) const
{
        if (!m_doc.save_file(path.c_str(), g_indent, pugi::format_default)) {
                throw std::runtime_error("stride::util::ptree::ToXmlFile> failed to write " + path);
        }
}

} // namespace util
} // namespace stride
