/*  Stream inserter for static_poly, with helpers.
 *  (C) Copyright Nick Matteo 2016.
 *  Use, modification and distribution are subject to the
 *  Boost Software License, Version 1.0. (See accompanying file
 *  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef NAM_STATIC_POLYNOMIAL_IO_HPP
#define NAM_STATIC_POLYNOMIAL_IO_HPP 

#include <ostream>
#include <cmath> //isnormal, fabs
#include <boost/range/algorithm/find_if.hpp>
#include <boost/range/algorithm/count_if.hpp>
#include <boost/math/special_functions/relative_difference.hpp>
#include "static_poly.hpp"

/** Forward declarations for ostream inserter helpers **/
namespace smath {
   template <typename T>
   struct complex;
}

namespace boost { namespace math {
   template <class T>
   class quaternion;

   template <class T>
   class octonion;
}}

/** Helpers for the stream inserter **/
namespace detail {
   struct xpow {
      int i;
   };

   inline std::ostream& operator << (std::ostream& os, xpow x) {
      if (x.i == 1)
         os << 'x';
      else if (x.i > 1)
         os << "x^" << x.i;
      return os;
   }

   template <typename T>
   std::enable_if_t<std::is_floating_point<T>::value, bool>
   is_zero(T n) {
       if (!std::isnormal(n))
           return true; // zero or denorm, or NaN or ∞
       return std::fabs(n) < 1e-11;
   }

   template <typename T>
   std::enable_if_t<!std::is_floating_point<T>::value, bool>
   is_zero(T n) {
       return n == T{0};
   }

   template <typename T>
   bool is_zero(smath::complex<T> ct) {
       return is_zero(ct.real) && is_zero(ct.imag.value);
   }

   template <typename T>
   struct notone {
       T val;
   };

   template <typename T>
   notone<T> ifnotone(T val) {
       return {val};
   }

   template <typename T>
   std::enable_if_t<std::is_floating_point<T>::value, bool>
   isone(T n) {
      using boost::math::relative_difference;
      return relative_difference(1.0, n) < 1e-11;
   }

   template <typename T>
   std::enable_if_t<!std::is_floating_point<T>::value, bool>
   isone(T n) {
      return n == T{1};
   }

   template <typename T>
   bool isone(smath::complex<T> n) {
      return isone(n.real) && is_zero(n.imag.value);
   }

   template <typename T>
   std::ostream& operator << (std::ostream& os, notone<T> io) {
      if (isone(-io.val))
         return os << '-';
      if (isone(io.val))
         return os;
      return os << io.val;
   }

   template <typename T>
   bool is_negative(T t) {
       return t < T{0} && !is_zero(t);
   }

   template <typename T>
   bool is_negative(smath::complex<T> ct) {
      return (ct.real < T{0} && !is_zero(ct.real)) ||
            (is_zero(ct.real) && ct.imag.value < T{0} && !is_zero(ct.imag.value));
   }

   template <typename T>
   bool is_negative(boost::math::quaternion<T> q) {
      // only if the first nonzero element is negative, and the majority of
      // nonzero terms are negative.
      using boost::range::find_if;
      using boost::range::count_if;
      T vals[] = {q.R_component_1(), q.R_component_2(),
                  q.R_component_3(), q.R_component_4()};
      auto nonzero = find_if(vals, [](T val){ return !is_zero(val); });
      if (nonzero == boost::end(vals)) return false;
      return *nonzero < T{0} &&
             count_if(vals, [](T val){ return val < T{0} && !is_zero(val); }) >=
             count_if(vals, [](T val){ return val > T{0} && !is_zero(val); });
   }

   template <typename T>
   bool is_negative(boost::math::octonion<T> q) {
      // only if the first nonzero element is negative, and the majority of
      // nonzero terms are negative.
      using boost::range::find_if;
      using boost::range::count_if;
      T vals[] = {q.R_component_1(), q.R_component_2(),
                  q.R_component_3(), q.R_component_4(),
                  q.R_component_5(), q.R_component_6(),
                  q.R_component_7(), q.R_component_8()};
      auto nonzero = find_if(vals, [](T val){ return !is_zero(val); });
      if (nonzero == boost::end(vals)) return false;
      return *nonzero < T{0} &&
             count_if(vals, [](T val){ return val < T{0} && !is_zero(val); }) >=
             count_if(vals, [](T val){ return val > T{0} && !is_zero(val); });
   }
} // close namespace detail

/** Output for smath::complex **/

namespace smath { // so it can be found
    template <typename T>
    std::ostream& operator << (std::ostream& os, const smath::complex<T>& cd) {
        using ::detail::is_zero;
        using ::detail::ifnotone;
        if (is_zero(cd.real)) {// real part is zero or denorm (or ∞, NaN)
            if (is_zero(cd.imag.value))
                return os << '0';
            return os << ifnotone(cd.imag.value) << 'i';
        }
        if (is_zero(cd.imag.value)) // imag part is zero or denorm (or ∞, NaN)
            return os << cd.real;
        os << '(' << cd.real;
        if (cd.imag.value < 0.)
            os << " - " << ifnotone(-cd.imag.value);
        else
            os << " + " << ifnotone(cd.imag.value);
        return os << "i)";
    }
}

template <class T, int N>
inline std::ostream& operator << (std::ostream& os, const static_poly<T, N>& poly) {
   using namespace detail;

   int i = poly.degree();
   if (i == -1)
      return os << '0';
   if (i == 0)
      return os << poly[0];

   os << ifnotone(poly[i]) << xpow{i};

   for (--i; i > 0; --i) {
      if (is_negative(poly[i]))
         os << " - " << ifnotone(-poly[i]) << xpow{i};
      else if (!is_zero(poly[i]))
         os << " + " << ifnotone(poly[i]) << xpow{i};
   }

   if (is_negative(poly[0]))
      os << " - " << -poly[0];
   else if (poly[0] != T{0})
      os << " + " << poly[0];
   return os;
}

#endif // NAM_STATIC_POLYNOMIAL_IO_HPP 

