# static-poly
Compile-time polynomials (C++14)

This allows creating and doing arithmetic with polynomials at compile time. For instance,

```C++
constexpr static_poly<int, 3> x2p1{1,0,1}; // x^2 + 1
constexpr static_poly<int, 2> xm1{-1,1}; // x - 1
constexpr auto prod = x2p1 * xm1; // x^3 - x^2 + x - 1
```

computes a polynomial product during compilation.

Currently, only clang++ seems capable of compiling this; I tried with clang-703.0.31. GCC 6.1 has an internal compiler error. GCC 5.3 says `'(&#'result_decl' not supported by dump_expr#<expression error>)`.
