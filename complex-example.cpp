#include <iostream>
#include <boost/math/constants/constants.hpp> // pi
#include <static_math/complex.h> // constexpr complex class
#include "static_poly_io.hpp"

using std::cout;

/*** Greatest common divisor ***/
constexpr int gcd(int m, int n) {
   while (n != 0) {
      m %= n;
      if (m == 0) return n;
      n %= m;
    }
    return m;
}

/*** The Nth cyclotomic polynomial, computed with complex numbers ***/
template <int N>
constexpr static_poly<smath::complex<double>, N+1> cyclotomic() {
    static_assert(N > 0, "No nonpositive numbers please!");
    using namespace boost::math::double_constants;

    static_poly<smath::complex<double>, N+1> cyc{1};

    for (int k = 1; k <= N; ++k) {
        if (gcd(k, N) == 1) {
            static_poly<smath::complex<double>, 2> x{-smath::polar(1., 2.*k*pi / N), 1};
            cyc = detail::mul(cyc, x);
        // alas, std::polar(1., 2.*k*pi / N) is not constexpr
        // furthermore, std::cos/std::sin are not constexpr (why!?)
        }
    }
    return cyc;
}

/*** Euler's totient function... ***/
template <int N>
constexpr int euler_totient() {
    return cyclotomic<N>().degree();
}

template <int N>
struct foo {
    constexpr int val() const {
        return N;
    }
};

int main() {
    constexpr auto phi1 = cyclotomic<1>();
    constexpr auto phi2 = cyclotomic<2>();
    constexpr auto phi3 = cyclotomic<3>();
    constexpr auto phi4 = cyclotomic<4>();
    constexpr auto phi5 = cyclotomic<5>();
    constexpr auto phi6 = cyclotomic<6>();
    constexpr auto phi7 = cyclotomic<7>();
    constexpr auto phi8 = cyclotomic<8>();
    constexpr auto phi9 = cyclotomic<9>();
    constexpr auto phi35 = cyclotomic<35>();

    cout << phi1 << '\n'
         << phi2 << '\n'
         << phi3 << '\n'
         << phi4 << '\n'
         << phi5 << '\n'
         << phi6 << '\n'
         << phi7 << '\n'
         << phi8 << '\n'
         << phi9 << "\n\n"
         << phi1*phi2 << '\n'
         << phi1*phi3 << '\n'
         << phi1*phi2*phi4 << '\n'
         << phi1*phi5 << '\n'
         << phi1*phi2*phi3*phi6 << '\n'
         << phi1*phi7 << '\n'
         << phi1*phi2*phi4*phi8 << '\n'
         << phi1*phi3*phi9 << "\n\n"
         << phi35 << '\n';

    /* Use euler's totient function as a template parameter! */
    foo<euler_totient<76>()> myfoo;
    cout << myfoo.val() << '\n';

    return 0;
}


