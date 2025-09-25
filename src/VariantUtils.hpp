/*
 * Authors: Gustavo V. Barroso
 * Created: 18/09/2025
 * Last modified: 25/09/2025
 *
 * Helpers to make Eigen compile safely
 */

#ifndef VARIANTUTILS_HPP
#define VARIANTUTILS_HPP

#include <iostream>
#include <variant>
#include <type_traits>
#include <cxxabi.h>
#include <Bpp/Exceptions.h>    // for bpp::Exception

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
// visitSameType(a,b,fn):
//   1) runtime-check that a.index()==b.index()
//   2) compile-time-check that the active types are identical
//   3) invoke fn(x,y) where x,y are the unwrapped values
// Throws bpp::Exception on mismatch.
//-----------------------------------------------------------------------------
template<class V1, class V2, class Fn>
void visitSameType(V1 &a, V2 const &b, Fn&& fn)
{
  if(a.index() != b.index())
    throw bpp::Exception("visitSameType: variant-index mismatch");

  #ifdef DEBUG
  auto demangle = [](const std::type_info& ti)
  {
    int status;
    char* demangled = abi::__cxa_demangle(ti.name(), nullptr, nullptr, &status);
    std::string result = (status == 0 && demangled) ? demangled : ti.name();
    free(demangled);
    return result;
  };
  #endif

  std::visit([&](auto &x, auto const &y)
  {
    using X = std::decay_t<decltype(x)>;
    using Y = std::decay_t<decltype(y)>;

    #ifdef DEBUG
    std::cout << "Inside visitSameType:\n";
    std::cout << "A type: " << demangle(typeid(X)) << ", address: " << static_cast<const void*>(&x) << "\n";
    std::cout << "B type: " << demangle(typeid(Y)) << ", address: " << static_cast<const void*>(&y) << "\n";
    #endif

   // NOTE: if constexpr(!std::is_same_v<decltype(x), decltype(y)>)
   // throws because it considers const, volatile and reference qualifiers (too restrictive)
   if constexpr(!std::is_same_v<X, Y>)
     throw bpp::Exception("visitSameType: scalar-type mismatch");

   else
    fn(x, y);

  }, a, b);
}

#endif // VARIANTUTILS_HPP
