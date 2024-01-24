/*
 * PCG Random Number Generation for C++
 *
 * Copyright 2014-2017 Melissa O'Neill <oneill@pcg-random.org>,
 *                     and the PCG Project contributors.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 * Licensed under the Apache License, Version 2.0 (provided in
 * LICENSE-APACHE.txt and at http://www.apache.org/licenses/LICENSE-2.0)
 * or under the MIT license (provided in LICENSE-MIT.txt and at
 * http://opensource.org/licenses/MIT), at your option. This file may not
 * be copied, modified, or distributed except according to those terms.
 *
 * Distributed on an "AS IS" BASIS, WITHOUT WARRANTY OF ANY KIND, either
 * express or implied.  See your chosen license for details.
 *
 * For additional information about the PCG random number generation scheme,
 * visit http://www.pcg-random.org/.
 */

/*
 * This file provides support code that is useful for random-number generation
 * but not specific to the PCG generation scheme, including:
 *      - 128-bit int support for platforms where it isn't available natively
 *      - bit twiddling operations
 *      - I/O of 128-bit and 8-bit integers
 *      - Handling the evilness of SeedSeq
 *      - Support for efficiently producing random numbers less than a given
 *        bound
 */

#ifndef PCG_EXTRAS_HPP_INCLUDED
#define PCG_EXTRAS_HPP_INCLUDED 1

namespace pcg_extras {

/*
 * We often need to represent a "number of bits".  When used normally, these
 * numbers are never greater than 128, so an unsigned char is plenty.
 * If you're using a nonstandard generator of a larger size, you can set
 * PCG_BITCOUNT_T to have it define it as a larger size.  (Some compilers
 * might produce faster code if you set it to an unsigned int.)
 */

#ifndef PCG_BITCOUNT_T
    typedef uint8_t bitcount_t;
#else
    typedef PCG_BITCOUNT_T bitcount_t;
#endif

/*
 * The C++ SeedSeq concept (modelled by seed_seq) can fill an array of
 * 32-bit integers with seed data, but sometimes we want to produce
 * larger or smaller integers.
 *
 * The following code handles this annoyance.
 *
 * uneven_copy will copy an array of 32-bit ints to an array of larger or
 * smaller ints (actually, the code is general it only needing forward
 * iterators).  The copy is identical to the one that would be performed if
 * we just did memcpy on a standard little-endian machine, but works
 * regardless of the endian of the machine (or the weirdness of the ints
 * involved).
 *
 * generate_to initializes an array of integers using a SeedSeq
 * object.  It is given the size as a static constant at compile time and
 * tries to avoid memory allocation.  If we're filling in 32-bit constants
 * we just do it directly.  If we need a separate buffer and it's small,
 * we allocate it on the stack.  Otherwise, we fall back to heap allocation.
 * Ugh.
 *
 * generate_one produces a single value of some integral type using a
 * SeedSeq object.
 */

 /* uneven_copy helper, case where destination ints are less than 32 bit. */

template<class SrcIter, class DestIter>
SrcIter uneven_copy_impl(
    SrcIter src_first, DestIter dest_first, DestIter dest_last,
    std::true_type)
{
    typedef typename std::iterator_traits<SrcIter>::value_type  src_t;
    typedef typename std::iterator_traits<DestIter>::value_type dest_t;

    constexpr bitcount_t SRC_SIZE  = sizeof(src_t);
    constexpr bitcount_t DEST_SIZE = sizeof(dest_t);
    constexpr bitcount_t DEST_BITS = DEST_SIZE * 8;
    constexpr bitcount_t SCALE     = SRC_SIZE / DEST_SIZE;

    size_t count = 0;
    src_t value = 0;

    while (dest_first != dest_last) {
        if ((count++ % SCALE) == 0)
            value = *src_first++;       // Get more bits
        else
            value >>= DEST_BITS;        // Move down bits

        *dest_first++ = dest_t(value);  // Truncates, ignores high bits.
    }
    return src_first;
}

 /* uneven_copy helper, case where destination ints are more than 32 bit. */

template<class SrcIter, class DestIter>
SrcIter uneven_copy_impl(
    SrcIter src_first, DestIter dest_first, DestIter dest_last,
    std::false_type)
{
    typedef typename std::iterator_traits<SrcIter>::value_type  src_t;
    typedef typename std::iterator_traits<DestIter>::value_type dest_t;

    constexpr auto SRC_SIZE  = sizeof(src_t);
    constexpr auto SRC_BITS  = SRC_SIZE * 8;
    constexpr auto DEST_SIZE = sizeof(dest_t);
    constexpr auto SCALE     = (DEST_SIZE+SRC_SIZE-1) / SRC_SIZE;

    while (dest_first != dest_last) {
        dest_t value(0UL);
        unsigned int shift = 0;

        for (size_t i = 0; i < SCALE; ++i) {
            value |= dest_t(*src_first++) << shift;
            shift += SRC_BITS;
        }

        *dest_first++ = value;
    }
    return src_first;
}

/* uneven_copy, call the right code for larger vs. smaller */

template<class SrcIter, class DestIter>
inline SrcIter uneven_copy(SrcIter src_first,
                           DestIter dest_first, DestIter dest_last)
{
    typedef typename std::iterator_traits<SrcIter>::value_type  src_t;
    typedef typename std::iterator_traits<DestIter>::value_type dest_t;

    constexpr bool DEST_IS_SMALLER = sizeof(dest_t) < sizeof(src_t);

    return uneven_copy_impl(src_first, dest_first, dest_last,
                            std::integral_constant<bool, DEST_IS_SMALLER>{});
}


template <size_t size, typename SeedSeq, typename DestIter>
void generate_to_impl(SeedSeq&& generator, DestIter dest,
                      std::false_type)
{
    typedef typename std::iterator_traits<DestIter>::value_type dest_t;
    constexpr auto DEST_SIZE = sizeof(dest_t);
    constexpr auto GEN_SIZE  = sizeof(uint32_t);

    constexpr bool GEN_IS_SMALLER = GEN_SIZE < DEST_SIZE;
    constexpr size_t FROM_ELEMS =
        GEN_IS_SMALLER
            ? size * ((DEST_SIZE+GEN_SIZE-1) / GEN_SIZE)
            : (size + (GEN_SIZE / DEST_SIZE) - 1)
                / ((GEN_SIZE / DEST_SIZE) + GEN_IS_SMALLER);
                        //  this odd code ^^^^^^^^^^^^^^^^^ is work-around for
                        //  a bug: http://llvm.org/bugs/show_bug.cgi?id=21287

    if (FROM_ELEMS <= 1024) {
        uint32_t buffer[FROM_ELEMS];
        generator.generate(buffer, buffer+FROM_ELEMS);
        uneven_copy(buffer, dest, dest+size);
    } else {
        uint32_t* buffer = static_cast<uint32_t*>(malloc(GEN_SIZE * FROM_ELEMS));
        generator.generate(buffer, buffer+FROM_ELEMS);
        uneven_copy(buffer, dest, dest+size);
        free(static_cast<void*>(buffer));
    }
}

template <size_t size, typename SeedSeq, typename DestIter>
inline void generate_to(SeedSeq&& generator, DestIter dest)
{
    typedef typename std::iterator_traits<DestIter>::value_type dest_t;
    constexpr bool IS_32BIT = sizeof(dest_t) == sizeof(uint32_t);

    generate_to_impl<size>(std::forward<SeedSeq>(generator), dest,
                           std::integral_constant<bool, IS_32BIT>{});
}

/* generate_one, produce a value of integral type using a SeedSeq
 * (optionally, we can have it produce more than one and pick which one
 * we want)
 */

template <typename UInt, size_t i = 0UL, size_t N = i+1UL, typename SeedSeq>
inline UInt generate_one(SeedSeq&& generator)
{
    UInt result[N];
    generate_to<N>(std::forward<SeedSeq>(generator), result);
    return result[i];
}


} // namespace pcg_extras

#endif // PCG_EXTRAS_HPP_INCLUDED
