/*
 * Authors: Gustavo V. Barroso
 * Created: 18/09/2025
 * Last modified: 25/09/2025
 *
 * Helpers to make Eigen compile safely
 *
 * This header preserves the original public intent:
 * - overloaded helper for combining lambdas
 * - always_false SFINAE helper
 * - demangle helper for debug diagnostics
 * - visitSameType: runtime index check, compile-time same-type guard, then call fn
 * - visitMatrixAndVector: runtime index check for matrix/vector variants and guarded dispatch
 *
 * Important: the visit functions use if constexpr guards so that user-provided
 * callables are only instantiated for valid type combinations.
 */

#ifndef VARIANTUTILS_HPP
#define VARIANTUTILS_HPP

#include <iostream>
#include <variant>
#include <type_traits>
#include <cxxabi.h>
#include <string>
#include <Bpp/Exceptions.h>
#include "Matrix.hpp" //  Matrix<Scalar> wrapper
#include "Vector.hpp" // Vector<Scalar> wrapper

//-----------------------------------------------------------------------------
// overloaded<Ts...>: helper to combine multiple lambdas into one callable
// usage:
//    std::visit(overloaded{
//      [&](A const& a){ ...  },
//      [&](B const& b){ ...  }
//    }, some_variant);
//-----------------------------------------------------------------------------
template<class... Ts>
struct overloaded : Ts...
{
  using Ts::operator()...;
};

// deduction guide
template<class... Ts>
overloaded(Ts...) -> overloaded<Ts...>;

//-----------------------------------------------------------------------------
// always_false<T>: SFINAE utility for static_assert in dependent contexts
//-----------------------------------------------------------------------------
template<typename T>
struct always_false : std::false_type
{};

//-----------------------------------------------------------------------------
// demangle helper (for DEBUG diagnostics)
//-----------------------------------------------------------------------------
inline std::string demangle_ti(const std::type_info& ti)
{
  int status = 0;
  char* demangled = abi::__cxa_demangle(ti.name(), nullptr, nullptr, &status);
  std::string result = (status == 0 && demangled) ? demangled : ti.name();
  free(demangled);
  return result;
}

//-----------------------------------------------------------------------------
// visitSameType(a,b,fn):
//   1) runtime-check that a.index()==b.index()
//   2) compile-time-check that the active unwrapped types are identical
//   3) invoke fn(x,y) where x,y are the unwrapped values
// Throws bpp::Exception on mismatch.
// Notes:
// - Does NOT require the two std::variant types themselves to be identical.
// - Uses if constexpr to avoid instantiating fn(x,y) for incompatible unwrapped types.
//-----------------------------------------------------------------------------
template<class V1, class V2, class Fn>
void visitSameType(V1 &a, V2 const &b, Fn&& fn)
{
  if(a.index() != b.index())
    throw bpp::Exception("visitSameType: variant-index mismatch");

  #ifdef DEBUG
  std::cout << "[visitSameType] variant index: " << a.index() << "\n";
  #endif

  std::visit([&](auto &x, auto const &y)
  {
    using X = std::decay_t<decltype(x)>;
    using Y = std::decay_t<decltype(y)>;

    #ifdef DEBUG
    std::cout << "Inside visitSameType:\n";
    std::cout << " A type: " << demangle_ti(typeid(X)) << ", addr: " << static_cast<const void*>(&x) << "\n";
    std::cout << " B type: " << demangle_ti(typeid(Y)) << ", addr: " << static_cast<const void*>(&y) << "\n";
    #endif

    if constexpr(std::is_same_v<X, Y>)
      fn(x, y); // instantiates and call fn only when unwrapped types match

    else
      throw bpp::Exception("visitSameType: scalar-type mismatch");

  }, a, b);
}

//-----------------------------------------------------------------------------
// visitMatrixAndVector(matrixVariant, vectorVariant, fn):
//   Convenience helper for dispatching matrix × vector calls when matrix and
//   vector are held in different variant types (e.g. MatrixVariant vs VectorVariant).
//   It checks runtime index compatibility (matrix type corresponds to vector type)
//   and then invokes fn(matrix, vector) guarded by a compile-time scalar-type check.
//   Throws bpp::Exception on mismatch.
//-----------------------------------------------------------------------------
template<class MatrixVariantT, class VectorVariantT, class Fn>
void visitMatrixAndVector(MatrixVariantT const &mvar, VectorVariantT &vvar, Fn&& fn)
{
  if(mvar.index() != vvar.index())
    throw bpp::Exception("visitMatrixAndVector: variant-index mismatch between matrix and vector");

  #ifdef DEBUG
  std::cout << "[visitMatrixAndVector] index: " << mvar.index() << "\n";
  #endif

  std::visit(overloaded
  {
    [&](const Matrix<double>& m, Vector<double>& v)
    {
      fn(m, v);
    },
    [&](const Matrix<mpfr::mpreal>& m, Vector<mpfr::mpreal>& v)
    {
      fn(m, v);
    },
    [&](auto const&, auto &) {
      throw bpp::Exception("visitMatrixAndVector: unsupported matrix/vector type combination");
    }
  }, mvar, vvar);
}

#endif // VARIANTUTILS_HPP
