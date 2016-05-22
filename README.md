# static-poly
Compile-time polynomials (C++14)

This allows creating and doing arithmetic with polynomials at compile time.
For instance,

```C++
constexpr static_poly<int, 3> x2p1{1,0,1}; // x^2 + 1
constexpr static_poly<int, 2> xm1{-1,1}; // x - 1
constexpr auto prod = x2p1 * xm1; // x^3 - x^2 + x - 1
```

computes a polynomial product during compilation.

Currently, only clang++ seems capable of compiling this; I tried with
clang-703.0.31.  `example.cpp` must be compiled with the option
`-fconstexpr-steps=4194304`.  GCC 5 and 6 have internal compiler errors;
see [PR71226](https://gcc.gnu.org/bugzilla/show_bug.cgi?id=71226).

The only requirements for the header `static_poly.hpp` itself are a C++14
standard libary and Boost (specifically Boost.MPL).  The `<<` operator for
stream output, in `static_poly_io.hpp`, also uses Boost.Range and Boost.Math.
These also suffice for the main example (`example.cpp`).

For the complex example, this [static_math library](https://github.com/kundor/static_math)
for compile-time mathematics support is required.  The standard library complex
class has constexpr constructors and `real()`, `imag()` accessors, but none of
the operators are constexpr!  This seems like something that could be changed
(it is [Library Working Group issue 2693](https://cplusplus.github.io/LWG/lwg-active.html#2693).)
