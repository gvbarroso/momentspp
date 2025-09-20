#pragma once
#include <variant>
#include <type_traits>
#include <utility>

namespace momentspp { namespace eigen_access {

// Single‐variant overloads
template<typename Variant, typename Fn>
void visit_eigen(Variant& var, Fn&& fn) {
  std::visit([&](auto& mat) -> void { fn(mat); }, var);
}
template<typename Variant, typename Fn>
void visit_eigen(const Variant& var, Fn&& fn) {
  std::visit([&](auto const& mat) -> void { fn(mat); }, var);
}

// Two‐variant overloads
template<typename V1, typename V2, typename Fn>
void visit_eigen(V1& v1, const V2& v2, Fn&& fn) {
  std::visit([&](auto& a, auto const& b) -> void {
    if constexpr(std::is_same_v<typename decltype(a)::Scalar,
                                typename decltype(b)::Scalar>)
      fn(a, b);
  }, v1, v2);
}
template<typename V1, typename V2, typename Fn>
void visit_eigen(V1& v1, V2& v2, Fn&& fn) {
  std::visit([&](auto& a, auto& b) -> void {
    if constexpr(std::is_same_v<typename decltype(a)::Scalar,
                                typename decltype(b)::Scalar>)
      fn(a, b);
  }, v1, v2);
}
template<typename V1, typename V2, typename Fn>
void visit_eigen(const V1& v1, const V2& v2, Fn&& fn) {
  std::visit([&](auto const& a, auto const& b) -> void {
    if constexpr(std::is_same_v<typename decltype(a)::Scalar,
                                typename decltype(b)::Scalar>)
      fn(a, b);
  }, v1, v2);
}

}} // namespace momentspp::eigen_access
