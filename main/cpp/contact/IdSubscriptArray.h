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
 *  Copyright 2019, Willem L, Kuylen E, Broeckhove J
 */

/**
 * @file
 * Container for the contact pools of various type (household, workplace, ...)..
 */

#pragma once

#include "contact/ContactType.h"

#include <array>
#include <bitset>
#include <initializer_list>
#include <iostream>
#include <stdexcept>

namespace stride {
namespace ContactType {

/**
 * An std::array modified to enable subscripting with the contact pool
 * type identifiers. The bounds checking "at" method, to protect against
 * a subscript that is static_cast<Id> from an erroneous integer value.
 *
 * We intentionally shadow the direct subscripting of the base array.
 *
 * Writing constructor is somewhat weird since std::array has no constructors,
 * only aggregate initialization. But we want some constructors for
 * shorthand expressions.
 *
 * @tparam T  base type of the array.
 */
template <class T>
class IdSubscriptArray : public std::array<T, NumOfTypes()>
{
public:
        /// What we will use most often and where we can have a default and
        /// initialize all array elements to the same value.
        /// e.g.    IdSubscriptArray<unsigned int> m(1U);
        explicit IdSubscriptArray(T t)
        {
                for (auto typ : IdList) {
                        this->operator[](typ) = t;
                }
        }

        explicit IdSubscriptArray()
        {
                for (auto typ : IdList) {
                        this->operator[](typ) = T();
                }
        }

        /// When we want to use an initializer list the elements is the
        /// (possibly empty) initializer list are applied to the first
        /// elements in order; any remaining elements are default initialized.
        ///  e.g.   IdSubscriptArray<unsigned int> mm {1U, 2U, 2U};
        IdSubscriptArray(std::initializer_list<T> l)
        {
                auto it = l.begin();
                for (auto typ : IdList) {
                        if (it != l.end()) {
                                this->operator[](typ) = *it;
                                ++it;
                        } else {
                                this->operator[](typ) = T();
                        }
                }
        }

        /// Initialize with an array of the right dimensions.
        /// Makes sure nothing is default constructed
        /// Delegate it to the array move constructor
        /// Destructive...
        explicit IdSubscriptArray(std::array<T, NumOfTypes()>&& l) : std::array<T, NumOfTypes()>(l) {}

        /// This actually works in itself but interferes annoyingly with the first
        /// constructor above and is for practical purpose redundant.
        ///
        /// Forward constructor arguments to base class.
        /// template <typename... Args>
        /// explicit IdSubscriptArray(Args&&... args) : std::array<T, NumOfTypes()>{{std::forward<Args>(args)...}}{}

        /// Subscripting with pool type id as argument.
        typename std::array<T, NumOfTypes()>::reference operator[](ContactType::Id id)
        {
                return this->std::template array<T, NumOfTypes()>::operator[](ContactType::ToSizeT(id));
        }

        /// Subscripting with pool type id as argument.
        typename std::array<T, NumOfTypes()>::const_reference operator[](ContactType::Id id) const
        {
                return this->std::template array<T, NumOfTypes()>::operator[](ContactType::ToSizeT(id));
        }

        /// Subscripting with pool type id as argument.
        typename std::array<T, NumOfTypes()>::reference at(ContactType::Id id)
        {
                if (ToSizeT(id) >= NumOfTypes()) {
                        throw std::out_of_range("IdSubscriptArray::at> Id out of range");
                }
                return this->std::template array<T, NumOfTypes()>::operator[](ContactType::ToSizeT(id));
        }

        /// Subscripting with pool type id as argument.
        typename std::array<T, NumOfTypes()>::const_reference at(ContactType::Id id) const
        {
                if (ToSizeT(id) >= NumOfTypes()) {
                        throw std::out_of_range("IdSubscriptArray::at> Id out of range");
                }
                return this->std::template array<T, NumOfTypes()>::operator[](ContactType::ToSizeT(id));
        }
};

/**
 * Specialization of IdSubscriptArray for booleans.
 */
template <>
class IdSubscriptArray<bool> : public std::bitset<NumOfTypes()>
{
public:
        /// What we will use most often and where we can have a default and
        /// initialize all array elements to the same value.
        /// e.g.    IdSubscriptArray<bool> m(true);
        explicit IdSubscriptArray(bool t = bool())
        {
                if (t)
                        this->set();
                else
                        this->reset();
        }

        /// When we want to use an initializer list the elements is the
        /// (possibly empty) initializer list are applied to the first
        /// elements in order; any remaining elements are default initialized.
        ///  e.g.   IdSubscriptArray<bool> mm {true, false, true};
        IdSubscriptArray(std::initializer_list<bool> l)
        {
                auto it = l.begin();
                for (auto typ : IdList) {
                        if (it != l.end()) {
                                this->operator[](typ) = *it;
                                ++it;
                        } else {
                                this->operator[](typ) = bool();
                        }
                }
        }

        /// Subscripting with pool type id as argument.
        typename std::bitset<NumOfTypes()>::reference operator[](ContactType::Id id)
        {
                return this->std::template bitset<NumOfTypes()>::operator[](ContactType::ToSizeT(id));
        }

        /// Subscripting with pool type id as argument.
        bool operator[](ContactType::Id id) const
        {
                return this->std::template bitset<NumOfTypes()>::operator[](ContactType::ToSizeT(id));
        }

        /// Subscripting with pool type id as argument.
        typename bitset<NumOfTypes()>::reference at(ContactType::Id id)
        {
                if (ToSizeT(id) >= NumOfTypes()) {
                        throw std::out_of_range("IdSubscriptArray<bool>::at> Id out of range");
                }
                return this->std::template bitset<NumOfTypes()>::operator[](ContactType::ToSizeT(id));
        }

        /// Subscripting with pool type id as argument.
        bool at(ContactType::Id id) const
        {
                if (ToSizeT(id) >= NumOfTypes()) {
                        throw std::out_of_range("IdSubscriptArray::at> Id out of range");
                }
                return this->std::template bitset<NumOfTypes()>::operator[](ContactType::ToSizeT(id));
        }
};

} // namespace ContactType
} // namespace stride
