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
#include <limits> // is_integer
#include <utility> // pair
#include <initializer_list>
#include <stdexcept>
#include "evaluate.hpp"

template <typename T>
class static_poly;

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
template <typename T, typename N>
typename disable_if_c<std::numeric_limits<T>::is_integer, void >::type
division_impl(static_poly<T> &q, static_poly<T> &u, const static_poly<T>& v, N n, N k)
{
    q[k] = u[n + k] / v[n];
    for (N j = n + k; j > k;)
    {
        j--;
        u[j] -= q[k] * v[j - k];
    }
}

template <class T, class N>
T integer_power(T t, N n)
{
   switch(n)
   {
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
   if(n & 1)
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
template <typename T, typename N>
typename enable_if_c<std::numeric_limits<T>::is_integer, void >::type
division_impl(static_poly<T> &q, static_poly<T> &u, const static_poly<T>& v, N n, N k)
{
    q[k] = u[n + k] * integer_power(v[n], k);
    for (N j = n + k; j > 0;)
    {
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
template <typename T>
std::pair< static_poly<T>, static_poly<T> >
division(static_poly<T> u, const static_poly<T>& v)
{
    assert(v.size() <= u.size());
    assert(v);
    assert(u);

    typedef typename static_poly<T>::size_type N;
    
    N const m = u.size() - 1, n = v.size() - 1;
    N k = m - n;
    static_poly<T> q;
    q.data().resize(m - n + 1);

    do
    {
        division_impl(q, u, v, n, k);
    }
    while (k-- != 0);
    u.data().resize(n);
    return std::make_pair(q, u);
}

} // namespace detail

/* Calculates a / b and a % b, returning the pair (quotient, remainder) together
 * because the same amount of computation yields both.
 * This function is not defined for division by zero: user beware.
 */
template <typename T>
std::pair< static_poly<T>, static_poly<T> >
quotient_remainder(const static_poly<T>& dividend, const static_poly<T>& divisor) {
    assert(divisor);
    if (dividend.size() < divisor.size())
        return std::make_pair(static_poly<T,0>(), dividend);
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
   explicit constexpr static_poly(const U (&data)[M]) : static_poly(&data, &data + M) {}

   template <class It>
   constexpr static_poly(It first, It last) {
      int i = 0;
      while (i < N && first != last)
         m_data[i++] = *first++;
      for (; i < N; ++i)
         m_data[i] = 0;
   }

   template <class U>
   explicit constexpr static_poly(const U& point) : m_data{point} {}

   // copy:
   static_poly(const static_poly& p) : static_poly(p.m_data) {}
   // call the array constructor; p has the same size as this.

   template <class U, int N1>
   explicit static_poly(const static_poly<U, N1>& p)
   : static_poly(&p.m_data, &p.m_data + N1) {} // call the pair-of-iterators constructor.
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
      if (std::numeric_limits<T>::is_integer) {
         for (T& i : m_data)
            i -= T(value * T(i / value));
      } else {
         for (T& i : m_data)
            i = 0; // note: std::fill, memset, etc. not constexpr
      }
      return *this;
   }

   static_poly& operator >>=(int n)
   {
       assert(n >= 0 && static_cast<unsigned>(n) <= m_data.size());
       m_data.erase(m_data.begin(), m_data.begin() + n);
       return *this;
   }

   static_poly& operator <<=(int n)
   {
       assert(n >= 0);
       m_data.insert(m_data.begin(), n, static_cast<T>(0));
       return *this;
   }

   explicit constexpr operator bool() const {
      return degree() >= 0;
   }
};


template <class T, int N, class U>
constexpr static_poly<T, std::max(N, 1)> operator + (static_poly<T, N> a, const U& b) {
   if (N) return a += b;
   return static_poly<T, 1>(b);
}

template <class T, int N, class U>
constexpr static_poly<T, std::max(N, 1)> operator - (static_poly<T, N> a, const U& b) {
   if(N) return a -= b;
   return static_poly<T, 1>(-b);
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
   return static_poly<T, 1>(a);
}

template <class U, class T, int N>
constexpr static_poly<T,std::max(N, 1)> operator - (const U& a, const static_poly<T, N>& b) {
   static_poly<T, std::max(N, 1)> result(a);
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
constexpr static_poly<T, N1 + N2> operator * (const static_poly<T, N1>& a, const static_poly<T, N2>& b) {
   static_poly<T, N1 + N2> prod;
   if (!a || !b) { // a or b is zero
      return prod;
   }
   for (int i = 0; i < N1; ++i)
      for (int j = 0; j < N2; ++j)
         prod[i+j] += a[i] * b[j];
   return prod;
}

template <class T, int N1, int N2>
constexpr static_poly<T, std::max(N1 - N2, 0)> operator / (const static_poly<T, N1>& a, const static_poly<T, N2>& b) {
   return quotient_remainder(a, b).first;
}

template <class T, int N1, int N2>
constexpr static_poly<T, std::min(N1, N2)> operator %(const static_poly<T, N1>& a, const static_poly<T, N2>& b) {
   return quotient_remainder(a, b).second;
}

template <class T, int N1, int N2>
constexpr bool operator == (const static_poly<T, N1> &a, const static_poly<T, N2> &b) {
   int n = a.degree();
   if (b.degree() != n) return false;
   for (int i = 0; i < n; ++i)
      if (a[i] != b[i]) return false;
   return true;
}

template <class T, int N1, int N2>
constexpr bool operator != (const static_poly<T, N1> &a, const static_poly<T, N2> &b) {
   return !(a == b);
}

template <class T, int N1, int N2>
constexpr bool operator < (const static_poly<T, N1> &a, const static_poly<T, N2> &b) {
    if (a.degree() != b.degree())
        return a.degree() < b.degree();
    return std::lexicographical_compare(a.data().rbegin(), a.data().rend(),
                                        b.data().rbegin(), b.data().rend());
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

template <class T, int N1>
constexpr auto operator >> (static_poly<T> a, int b) -> static_poly<T, std:::max(N - b, 0)> {

    a >>= b;
    return a;
}

template <class T>
constexpr static_poly<T> operator << (static_poly<T> a, int b) {
    a <<= b;
    return a;
}

// Unary minus (negate).
template <class T, int N>
constexpr static_poly<T, N> operator - (static_poly<T, N> a) {
   for (T& i : a.m_data)
      i *= -1;
   return a;
}

template <class T>
constexpr bool odd(static_poly<T> const &a) {
    return a.size() > 0 && a[0] != static_cast<T>(0);
}

template <class T>
constexpr bool even(static_poly<T> const &a) {
    return !odd(a);
}

template <class T>
constexpr static_poly<T> pow(static_poly<T> base, int exp) {
    if (exp < 0)
        throw(std::domain_error("Negative powers are not supported for polynomials.");
        // if this is called as constexpr, trying to throw is a compile-time error: good
    static_poly<T> result(T(1));
    if (exp & 1)
        result = base;
    /* "Exponentiation by squaring" */
    while (exp >>= 1) {
        base *= base;
        if (exp & 1)
            result *= base;
    }
    return result;
}

template <class charT, class traits, class T>
inline std::basic_ostream<charT, traits>& operator << (std::basic_ostream<charT, traits>& os, const static_poly<T>& poly) {
   os << "{ ";
   for(unsigned i = 0; i < poly.size(); ++i)
   {
      if(i) os << ", ";
      os << poly[i];
   }
   os << " }";
   return os;
}

#endif // NAM_STATIC_POLYNOMIAL_HPP



