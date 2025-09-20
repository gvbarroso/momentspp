#ifndef VISIT_SAME_SCALAR_HPP
#define VISIT_SAME_SCALAR_HPP

#include <variant>
#include <type_traits>
#include "MatrixEngine.hpp"
#include <bpp/Exceptions.h>   // for bpp::Exception

/**
 *  Ensures we only call op(A, v) when A and v share the same Scalar.
 *  Throws bpp::Exception on mismatch.
 */
template <class MatVar, class VecVar, class Op>
void visit_same_scalar(MatVar&& matVar, VecVar&& vecVar, Op&& op) {
  std::visit([&](auto&& A) {
    using MScalar = typename std::decay_t<decltype(A)>::Scalar;
    std::visit([&](auto&& v) {
      using VScalar = typename std::decay_t<decltype(v)>::Scalar;
      if constexpr (std::is_same_v<MScalar, VScalar>) {
        op(A, v);
      } else {
        throw bpp::Exception("Scalar mismatch between matrix and vector variant");
      }
    }, std::forward<VecVar>(vecVar));
  }, std::forward<MatVar>(matVar));
}

#endif // VISIT_SAME_SCALAR_HPP
