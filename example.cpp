#include <iostream>
#include "static_poly.hpp"

using std::cout;

int main() {
    constexpr static_poly<int, 3> x2p1{1,0,1}; // x^2 + 1
    constexpr static_poly<int, 2> xm1{-1,1}; // x - 1
    constexpr auto prod = x2p1 * xm1; // x^3 - x^2 + x - 1
    cout << prod << '\n';
    return 0;
}

