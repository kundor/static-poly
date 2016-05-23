#include <iostream>
#include "static_poly.hpp"

using std::cout;

void foo(static_poly<int, 4> bar) {
    cout << "Foo: " << bar*2 << '\n';
}

constexpr int qux(static_poly<int, 4> bar) {
    return bar.degree();
}


int main() {
    constexpr static_poly<int, 3> x2p1{1,0,1}; // x^2 + 1
    constexpr static_poly<int, 2> xm1{-1,1}; // x - 1
    constexpr auto prod = x2p1 * xm1; // x^3 - x^2 + x - 1
    cout << prod << '\n';

    foo(prod); // succeeds, but if foo took a static_poly<int, 5>, would fail:
               // copying to other sizes is explicit

    constexpr static_poly<int, 6> cp{prod}; // succeeds: explicitly change size
    cout << cp << '\n';
    constexpr static_poly<int, 3> cp2{prod}; // succeeds: explicitly change size
    cout << cp2 << '\n';

    constexpr int i = qux(prod); // The default copy constructor is constexpr
    cout << "Degree" << i << "\n\n";

    constexpr static_poly<int, 2> x{0,1}; // x
    constexpr auto bigprod = (power<5>(x) - 3*power<4>(x) + 2*power<3>(x) - 7*x*x + 4*x + 1)
                           * (x*x*x - 6*x*x + x - 1)
                           * (x*x - 3);
    cout << bigprod << "\n"
         "Cubed:  " << power<3>(bigprod) << "\n\n";

    /* Cyclotomic polynomials */
    constexpr auto phi1 = x - 1;
    constexpr auto phi2 = x + 1;
    constexpr auto phi3 = x*x + x + 1;
    constexpr auto phi4 = x*x + 1;
    constexpr auto phi5 = (power<5>(x) - 1)/phi1; // x^4 + x^3 + x^2 + x + 1
    constexpr auto phi6 = x*x - x + 1;
    constexpr auto phi7 = (power<7>(x) - 1)/phi1; // x^6 + x^5 + x^4 + x^3 + x^2 + x + 1

    cout << "ϕ5: " << phi5 << '\n';
    cout << "ϕ7: " << phi7 << '\n';
    cout << "ϕ14: " << (power<14>(x) - 1)/phi7/phi2/phi1 << "\n\n";

    cout << phi1*phi2 << '\n'; // x^2 - 1
    cout << phi1*phi3 << '\n'; // x^3 - 1
    cout << phi1*phi2*phi4 << '\n'; // x^4 - 1
    cout << phi1*phi2*phi3*phi6 << "\n\n"; // x^6 - 1

    constexpr auto phi15 = (power<15>(x) - 1)/phi5/phi3/phi1;
    constexpr auto phi21 = (power<21>(x) - 1)/phi7/phi3/phi1;
    constexpr auto phi35 = (power<35>(x) - 1)/phi7/phi5/phi1;
    constexpr auto phi105 = (power<105>(x) - 1)/phi35/phi21/phi15/phi7/phi5/phi3/phi1;
    // in clang by default, "constexpr evaluation hit maximum step limit"
    // compile with -fconstexpr-steps=4194304
    cout << "ϕ105: " << phi105 << '\n'; // The first example with a coefficient other than ±1
    
    return 0;
}

