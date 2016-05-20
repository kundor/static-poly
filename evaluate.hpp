// Modified from boost/math/tools/rational.hpp (as of Boost 1.61)
//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef NAM_EVALUATE_HPP
#define NAM_EVALUATE_HPP

#include <boost/mpl/int.hpp>
#include "polynomial_horner2_20.hpp"

//
// Forward declaration to keep two phase lookup happy:
//
template <class T, class U>
constexpr U evaluate_polynomial(const T* poly, U const& z, std::size_t count);

namespace detail {
    // Fallback if none of the unrolled ones (up to size 20) suffice
    template <class T, class V, class Tag>
    constexpr V evaluate_polynomial_c_imp(const T* a, const V& val, const Tag*) {
       return evaluate_polynomial(a, val, Tag::value);
    }
} // namespace detail

//
// Polynomial evaluation with runtime size.
// This requires a for-loop which may be more expensive than
// the loop expanded versions above:
//
template <class T, class U>
constexpr U evaluate_polynomial(const T* poly, U const& z, std::size_t count) {
   assert(count > 0);
   U sum = static_cast<U>(poly[count - 1]);
   for(int i = static_cast<int>(count) - 2; i >= 0; --i) {
      sum *= z;
      sum += static_cast<U>(poly[i]);
   }
   return sum;
}

//
// Compile time sized polynomials, just inline forwarders to the
// implementations above:
//
template <std::size_t N, class T, class V>
constexpr V evaluate_polynomial(const T(&a)[N], const V& val) {
   typedef mpl::int_<N> tag_type;
   return detail::evaluate_polynomial_c_imp(static_cast<const T*>(a), val, static_cast<tag_type const*>(0));
}

#endif // NAM_EVALUATE_HPP




