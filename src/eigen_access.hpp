#pragma once

// momentspp/eigen_access.hpp
// Small utilities to normalize access to underlying Eigen objects.
// Provides: has_eigen_getter, eigen_ref, eigen_cref, visit_eigen overloads.

#include <type_traits>
#include <utility>
#include <variant>

namespace momentspp::eigen_access {

// detect .eigen() member
template<typename T, typename = void>
struct has_eigen_getter : std::false_type {};

template<typename T>
struct has_eigen_getter<T, std::void_t<decltype(std::declval<T>().eigen())>> : std::true_type {};

// eigen_ref for non-const objects: returns reference to underlying Eigen object
template<typename T>
inline auto eigen_ref(T& obj)
  -> std::enable_if_t<has_eigen_getter<T>::value, decltype(obj.eigen())>
{
  return obj.eigen();
}

template<typename T>
inline auto eigen_ref(T& obj)
  -> std::enable_if_t<!has_eigen_getter<T>::value, T&>
{
  return obj;
}

// eigen_cref for const objects: returns const reference to underlying Eigen object
template<typename T>
inline auto eigen_cref(const T& obj)
  -> std::enable_if_t<has_eigen_getter<T>::value, decltype(obj.eigen())>
{
  return obj.eigen();
}

template<typename T>
inline auto eigen_cref(const T& obj)
  -> std::enable_if_t<!has_eigen_getter<T>::value, const T&>
{
  return obj;
}

// visit_eigen single-variant (non-const)
template<typename Variant, typename Fn>
inline decltype(auto) visit_eigen(Variant& v, Fn&& fn)
{
  return std::visit([&](auto& alt) -> decltype(auto)
  {
    return std::forward<Fn>(fn)(eigen_ref(alt));
  }, v);
}

// visit_eigen single-variant (const)
template<typename Variant, typename Fn>
inline decltype(auto) visit_eigen(const Variant& v, Fn&& fn)
{
  return std::visit([&](const auto& alt) -> decltype(auto)
  {
    return std::forward<Fn>(fn)(eigen_cref(alt));
  }, v);
}

// visit_eigen two-variant overload (non-const)
template<typename V1, typename V2, typename Fn>
inline decltype(auto) visit_eigen(V1& v1, V2& v2, Fn&& fn)
{
  return std::visit([&](auto& a, auto& b) -> decltype(auto)
  {
    return std::forward<Fn>(fn)(eigen_ref(a), eigen_ref(b));
  }, v1, v2);
}

// visit_eigen two-variant overload (mix const/non-const as needed)
template<typename V1, typename V2, typename Fn>
inline decltype(auto) visit_eigen(const V1& v1, V2& v2, Fn&& fn)
{
  return std::visit([&](const auto& a, auto& b) -> decltype(auto)
  {
    return std::forward<Fn>(fn)(eigen_cref(a), eigen_ref(b));
  }, v1, v2);
}

} // namespace momentspp::eigen_access
