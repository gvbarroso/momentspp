/*
 * Authors: Gustavo V. Barroso
 * Created: 10/09/2025
 * Last modified: 10/09/2025
 *
 */

// NOTE this header is not used at the moment, ie, not included in any file

#ifndef EIGEN_MPREAL_TRAITS_HPP
#define EIGEN_MPREAL_TRAITS_HPP

#include <mpreal.h>            // MPFR C++ wrapper
#include <eigen3/Eigen/Core>   // Core Eigen functionality
#include <eigen3/Eigen/Sparse> // Sparse matrix support

// this helps Eigen use mpreal as a scalar
namespace Eigen
{
template <>
struct NumTraits<mpfr::mpreal> : GenericNumTraits<mpfr::mpreal>
{
  typedef mpfr::mpreal Real;
  typedef mpfr::mpreal NonInteger;
  typedef mpfr::mpreal Nested;

  enum
  {
    IsComplex = 0,
    IsInteger = 0,
    IsSigned = 1,
    RequireInitialization = 1,
    ReadCost = 1,
    AddCost = 2,
    MulCost = 3
  };

  static inline Real epsilon()
  {
    return mpfr::mpreal::get_epsilon();
  }
  static inline Real dummy_precision()
  {
    return mpfr::mpreal::get_epsilon();
  }
  static inline int digits10()
  {
    return mpfr::mpreal::get_default_prec();
  }
};
} // namespace Eigen

namespace std
{
inline mpfr::mpreal abs(const mpfr::mpreal& x)
{
  return mpfr::abs(x);
}
inline mpfr::mpreal sqrt(const mpfr::mpreal& x)
{
  return mpfr::sqrt(x);
}
inline mpfr::mpreal exp(const mpfr::mpreal& x)
{
  return mpfr::exp(x);
}
inline mpfr::mpreal log(const mpfr::mpreal& x)
{
  return mpfr::log(x);
}
inline mpfr::mpreal pow(const mpfr::mpreal& x, const mpfr::mpreal& y)
{
  return mpfr::pow(x, y);
}
inline mpfr::mpreal sin(const mpfr::mpreal& x)
{
  return mpfr::sin(x);
}
inline mpfr::mpreal cos(const mpfr::mpreal& x)
{
  return mpfr::cos(x);
}
} // namespace std

#endif
