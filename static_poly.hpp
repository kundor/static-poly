//  (C) Copyright Nick Matteo 2016.
//  Adapted from boost/math/tools/polynomial.hpp (from Boost 1.61) containing these notices:

//  (C) Copyright John Maddock 2006.
//  (C) Copyright Jeremy William Murphy 2015.


//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef NAM_STATIC_POLYNOMIAL_HPP
#define NAM_STATIC_POLYNOMIAL_HPP

#include <cassert>
#include <ostream>
#include <algorithm> // minmax
#include <type_traits> // enable_if, is_integral
#include <utility> // pair
#include <initializer_list>
#include <boost/range/algorithm/find_if.hpp>
#include <boost/range/algorithm/count_if.hpp>
#include "evaluate.hpp"

template <typename T, int N>
struct static_poly;

namespace detail {

/**
* Knuth, The Art of Computer Programming: Volume 2, Third edition, 1998
* Chapter 4.6.1, Algorithm D: Division of polynomials over a field.
*
* @tparam  T   Coefficient type, must be not be an integer.
*
* Template-parameter T actually must be a field but we don't currently have that
* subtlety of distinction.
*/
template <typename T, int N1, int N2, int N3>
std::enable_if_t<!std::is_integral<T>::value> /*void*/
constexpr division_impl(static_poly<T, N3> &q, static_poly<T, N1> &u, const static_poly<T, N2>& v, int n, int k) {
    q[k] = u[n + k] / v[n];
    for (int j = n + k; j > k;) {
        j--;
        u[j] -= q[k] * v[j - k];
    }
}

template <class T>
constexpr T integer_power(T t, int n) {
   switch(n) {
   case 0:
      return static_cast<T>(1u);
   case 1:
      return t;
   case 2:
      return t * t;
   case 3:
      return t * t * t;
   }
   T result = integer_power(t, n / 2);
   result *= result;
   if (n & 1)
      result *= t;
   return result;
}


/**
* Knuth, The Art of Computer Programming: Volume 2, Third edition, 1998
* Chapter 4.6.1, Algorithm R: Pseudo-division of polynomials.
*
* @tparam  T   Coefficient type, must be an integer.
*
* Template-parameter T actually must be a unique factorization domain but we
* don't currently have that subtlety of distinction.
*/
template <typename T, int N1, int N2, int N3>
std::enable_if_t<std::is_integral<T>::value> /*void*/
constexpr division_impl(static_poly<T, N3> &q, static_poly<T, N1> &u, const static_poly<T, N2>& v, int n, int k) {
   q[k] = u[n + k] * integer_power(v[n], k);
   for (int j = n + k; j > 0;) {
      j--;
      u[j] = v[n] * u[j] - (j < k ? T(0) : u[n + k] * v[j - k]);
   }
}


/**
 * Knuth, The Art of Computer Programming: Volume 2, Third edition, 1998
 * Chapter 4.6.1, Algorithm D and R: Main loop.
 *
 * @param   u   Dividend.
 * @param   v   Divisor.
 */
template <typename T, int N1, int N2>
std::pair< static_poly<T, std::max(N1 - N2 + 1, 1)>, static_poly<T, std::min(N1, N2)> >
constexpr division(static_poly<T, N1> u, const static_poly<T, N2>& v) {
   assert(v.size() <= u.size());
   assert(v);
   assert(u);

   const int m = u.degree(), n = v.degree();
   int k = m - n;
   static_poly<T, std::max(N1 - N2 + 1, 1)> q;

   do division_impl(q, u, v, n, k);
   while (k--); // stops when k was already 0

   return std::make_pair(q, static_poly<T, std::min(N1, N2)>(u));
}

} // namespace detail

/* Calculates a / b and a % b, returning the pair (quotient, remainder) together
 * because the same amount of computation yields both.
 * This function is not defined for division by zero: user beware.
 */
template <typename T, int N1, int N2>
std::pair< static_poly<T, std::max(N1 - N2 + 1, 1)>, static_poly<T, std::min(N1, N2)> >
constexpr quotient_remainder(const static_poly<T, N1>& dividend, const static_poly<T, N2>& divisor) {
   assert(divisor);
   constexpr int sz = std::max(N1 - N2 + 1, 1);
   if (dividend.degree() < divisor.degree())
      return std::make_pair(static_poly<T, sz>(), static_poly<T, std::min(N1, N2)>(dividend));
   return detail::division(dividend, divisor);
}


template <class T, int N>
struct static_poly {
   T m_data[N]; //constexpr std::array isn't modifiable in C++14 (P0107R0)

   // typedefs:
   typedef T value_type;
   typedef int size_type;

   // construct:
   constexpr static_poly() : m_data{} {}

   template <class U, int M>
   explicit constexpr static_poly(const U (&data)[M]) : static_poly(data, data + M) {}

   template <class It>
   constexpr static_poly(It first, It last) : m_data{} {
      int i = 0;
      while (i < N && first != last)
         m_data[i++] = *first++;
   }

   template <class U>
   explicit constexpr static_poly(const U& point) : m_data{point} {}

   // copy defaulted

   template <class U, int N1>
   explicit constexpr static_poly(const static_poly<U, N1>& p)
   : static_poly(p.m_data, p.m_data + N1) {} // call the pair-of-iterators constructor.
   // if N1 > N, we lose the high terms
   
   constexpr static_poly(std::initializer_list<T> l)
   : static_poly(std::begin(l), std::end(l)) {}

   // access:
   constexpr size_type size() const {
      return N;
   }

   constexpr size_type degree() const {
      for (int i = N-1; i >= 0; --i)
         if (m_data[i])
            return i;
      return -1; // zero polynomial: should really be -∞, or undefined.
   }

   constexpr T& operator[] (size_type i) {
      return m_data[i];
   }
   
   constexpr const T& operator[] (size_type i) const {
      return m_data[i];
   }
   
   constexpr T operator() (T z) const {
      return evaluate_polynomial(m_data, z);
   }
   
   constexpr const std::pair<const T*, const T*> data() const {
      // return a pair of iterators, suitable for use with Boost.Range
      return std::make_pair(&m_data, &m_data + N);
   }

   constexpr std::pair<T*, T*> & data() {
      return std::make_pair(&m_data, &m_data + N);
   }

   // operators:
   template <class U>
   constexpr static_poly& operator +=(const U& value) {
      static_assert(N, "Cannot modify zero polynomial");
      m_data[0] += value;
      return *this;
   }

   template <class U>
   constexpr static_poly& operator -=(const U& value) {
      static_assert(N, "Cannot modify zero polynomial");
      m_data[0] -= value;
      return *this;
   }

   template <class U>
   constexpr static_poly& operator *=(const U& value) {
      for (T& i : m_data)
         i *= value;
      return *this;
   }

   template <class U>
   constexpr static_poly& operator /=(const U& value) {
      for (T& i : m_data)
         i /= value;
      return *this;
   }

   template <class U>
   constexpr static_poly& operator %=(const U& value) {
      // In the case that T is integral, this preserves the semantics
      // p == r*(p/r) + (p % r), for polynomial<T> p and U r.
      if (std::is_integral<T>::value) {
         for (T& i : m_data)
            i -= T(value * T(i / value));
      } else {
         for (T& i : m_data)
            i = 0; // note: std::fill, memset, etc. not constexpr
      }
      return *this;
   }

   explicit constexpr operator bool() const {
      return degree() >= 0;
   }
};


template <class T, int N, class U>
constexpr static_poly<T, std::max(N, 1)> operator + (static_poly<T, N> a, const U& b) {
   if (N) return a += b;
   return static_poly<T, std::max(N, 1)>(b); // N is 0, so the max is 1
}

template <class T, int N, class U>
constexpr static_poly<T, std::max(N, 1)> operator - (static_poly<T, N> a, const U& b) {
   if(N) return a -= b;
   return static_poly<T, std::max(N, 1)>(-b);
}

template <class T, int N, class U>
constexpr static_poly<T, N> operator * (static_poly<T, N> a, const U& b) {
   return a *= b;
}

template <class T, int N, class U>
constexpr static_poly<T, N> operator / (static_poly<T, N> a, const U& b) {
   return a /= b;
}

template <class T, int N, class U>
constexpr static_poly<T, N> operator % (static_poly<T, N> a, const U& b) {
   return a %= b;
}

template <class U, class T, int N>
constexpr static_poly<T, std::max(N, 1)> operator + (const U& a, static_poly<T, N> b) {
   if (N) return b += a;
   return static_poly<T, std::max(N, 1)>(a);
}

template <class U, class T, int N>
constexpr static_poly<T,std::max(N, 1)> operator - (const U& a, const static_poly<T, N>& b) {
   static_poly<T, std::max(N, std::max(N, 1))> result(a);
   return result - b;
}

template <class U, class T, int N>
constexpr static_poly<T, N> operator * (const U& a, static_poly<T, N> b) {
   return b *= a;
}

template <class T, int N1, int N2>
constexpr static_poly<T, std::max(N1, N2)> operator + (const static_poly<T, N1>& a, const static_poly<T, N2>& b) {
   static_poly<T, std::max(N1, N2)> sum(a); // copies a's coefficients; if N2>N1, extends with 0
   for (int i = 0; i < N2; ++i)
      sum[i] += b[i];
   return sum;
}

template <class T, int N1, int N2>
constexpr static_poly<T, std::max(N1, N2)> operator - (const static_poly<T, N1>& a, const static_poly<T, N2>& b) {
   static_poly<T, std::max(N1, N2)> diff(a); // copies a's coefficients; if N2>N1, extends with 0
   for (int i = 0; i < N2; ++i)
      diff[i] -= b[i];
   return diff;
}

template <class T, int N1, int N2>
constexpr static_poly<T, N1 + N2 - 1> operator * (const static_poly<T, N1>& a, const static_poly<T, N2>& b) {
   static_poly<T, N1 + N2 - 1> prod;
   if (!a || !b) { // a or b is zero
      return prod;
   }
   for (int i = 0; i < N1; ++i)
      for (int j = 0; j < N2; ++j)
         prod[i+j] += a[i] * b[j];
   return prod;
}

namespace detail {
   // special case of multiplication: two polys into another of the same size as the first,
   // assuming that there's enough "headroom" (tail of zero coefficients) so the result fits
   template <class T, int N, int N2>
   constexpr static_poly<T, N> mul(const static_poly<T, N>& a, const static_poly<T, N2>& b) {
      static_poly<T, N> prod;
      if (!a || !b) { // a or b is zero
         return prod;
      }
      for (int i = 0; i < N; ++i)
         for (int j = 0; j < std::min(N - i, N2); ++j)
            prod[i+j] += a[i] * b[j];
      return prod;
   }
}

template <class T, int N1, int N2>
constexpr static_poly<T, std::max(N1 - N2 + 1, 1)> operator / (const static_poly<T, N1>& a, const static_poly<T, N2>& b) {
   return quotient_remainder(a, b).first;
}

template <class T, int N1, int N2>
constexpr static_poly<T, std::min(N1, N2)> operator % (const static_poly<T, N1>& a, const static_poly<T, N2>& b) {
   return quotient_remainder(a, b).second;
}

template <class T, int N1, int N2>
constexpr bool operator == (const static_poly<T, N1> &a, const static_poly<T, N2> &b) {
   int n = a.degree();
   if (b.degree() != n) return false;
   for (; n >= 0; --n)
      if (a[n] != b[n]) return false;
   return true;
}

template <class T, int N1, int N2>
constexpr bool operator != (const static_poly<T, N1> &a, const static_poly<T, N2> &b) {
   return !(a == b);
}

template <class T, int N1, int N2>
constexpr bool operator < (const static_poly<T, N1> &a, const static_poly<T, N2> &b) {
   int k = a.degree();
   if (b.degree() != k)
        return k < b.degree();
   for (; k >= 0; --k) {
      if (a[k] != b[k])
         return a[k] < b[k];
   }
   return false; // equal
}

template <class T, int N1, int N2>
constexpr bool operator <= (const static_poly<T, N1> &a, const static_poly<T, N2> &b) {
   return a < b || a == b;
}
   
template <class T, int N1, int N2>
constexpr bool operator >= (const static_poly<T, N1> &a, const static_poly<T, N2> &b) {
   return !(a < b);
}
 
template <class T, int N1, int N2>
constexpr bool operator > (const static_poly<T, N1> &a, const static_poly<T, N2> &b) {
   return !(a <= b);
}

// Unary minus (negate).
template <class T, int N>
constexpr static_poly<T, N> operator - (static_poly<T, N> a) {
   for (T& i : a.m_data)
      i *= -1;
   return a;
}

template <int exp, class T, int N>
constexpr static_poly<T, N*exp> power(const static_poly<T, N>& b) {
   static_assert(exp >= 0, "Negative power not supported");
   static_poly<T, N*exp> result{T{1}};
   static_poly<T, N*exp> base{b};
   int ex = exp;
   if (exp & 1)
        result = base;
    /* "Exponentiation by squaring" */
    while (ex >>= 1) {
        base = detail::mul(base, base);
        if (ex & 1)
            result = detail::mul(result, base);
    }
    return result;
}
/* A pow with an argument for exp would be nice--
 * but how big must the result be? We can't determine the return type
 * at compile time. */

/* Forward declarations for ostream inserter helpers */
namespace std {
#ifdef _LIBCPP_VERSION
// in libc++
  inline namespace __1 {
#endif

   template <typename T>
   class complex;
#ifdef _LIBCPP_VERSION
  }
#endif
}

namespace boost { namespace math {
   template <class T>
   class quaternion;

   template <class T>
   class octonion;
}}

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
   bool is_negative(T t) {
       return t < T{0};
   }

   template <typename T>
   bool is_negative(std::complex<T> ct) {
      return ct.real() < T{0} ||
            (ct.real() == T{0} && ct.imag() < T{0});
   }

   template <typename T>
   bool is_negative(boost::math::quaternion<T> q) {
      // only if the first nonzero element is negative, and the majority of
      // nonzero terms are negative.
      using boost::range::find_if;
      using boost::range::count_if;
      T vals[] = {q.R_component_1(), q.R_component_2(),
                  q.R_component_3(), q.R_component_4()};
      auto nonzero = find_if(vals, [](T val){ return val != T{0}; });
      if (nonzero == boost::end(vals)) return false;
      return *nonzero < T{0} &&
             count_if(vals, [](T val){ return val < T{0}; }) >=
             count_if(vals, [](T val){ return val > T{0}; });
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
      auto nonzero = find_if(vals, [](T val){ return val != T{0}; });
      if (nonzero == boost::end(vals)) return false;
      return *nonzero < T{0} &&
             count_if(vals, [](T val){ return val < T{0}; }) >=
             count_if(vals, [](T val){ return val > T{0}; });
   }
}

template <class T, int N>
inline std::ostream& operator << (std::ostream& os, const static_poly<T, N>& poly) {
   using detail::xpow;
   using detail::is_negative;
   int i = poly.degree();
   if (i == -1)
      return os << '0';
   if (i == 0)
      return os << poly[0];

   if (poly[i] == T{-1})
      os << '-';
   else if (poly[i] != T{1})
      os << poly[i];
   os << xpow{i};

   for (--i; i > 0; --i) {
      if (poly[i] == T{1})
         os << " + " << xpow{i};
      else if (poly[i] == T{-1})
         os << " - " << xpow{i};
      else if (is_negative(poly[i]))
         os << " - " << -poly[i] << xpow{i};
      else if (poly[i] != T{0})
         os << " + " << poly[i] << xpow{i};
   }

   if (is_negative(poly[0]))
      os << " - " << -poly[0];
   else if (poly[0] != T{0})
      os << " + " << poly[0];
   return os;
}
#endif // NAM_STATIC_POLYNOMIAL_HPP



