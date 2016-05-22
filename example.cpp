#include <iostream>
#include "static_poly.hpp"

using std::cout;

int main() {
    static_poly<int, 3> x2p1{1,0,1}; // x^2 + 1
    cout << x2p1 << '\n';
    return 0;
}

