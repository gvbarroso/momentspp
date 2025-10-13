/*
 * Authors: Gustavo V. Barroso
 * Created: 08/09/2025
 * Last modified: 13/10/2025
 *
 */

#ifndef _VECINTERFACE_H_
#define _VECINTERFACE_H_

#pragma once

#include <Bpp/Exceptions.h>

#include <eigen3/Eigen/Core>
#include <eigen3/unsupported/Eigen/MPRealSupport>
#include <mpreal.h>
#include <memory>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdexcept>

template <typename T>
class Vector
{
private:
  Eigen::Matrix<T, Eigen::Dynamic, 1> vec_;

public:
  using Scalar = T;
  using EigenVector = Eigen::Matrix<T, Eigen::Dynamic, 1>;

  Vector():
  vec_()
  { }

  explicit Vector(size_t size):
  vec_(size)
  {
    setZero();
  }

  Vector(const Vector& other):
  vec_(other.vec_)
  { }

  Vector(Eigen::Matrix<T, Eigen::Dynamic, 1>&& other):
  vec_(std::move(other))
  { }

  Vector<T>& operator=(Vector<T>&& other) noexcept
  {
    vec_ = std::move(other.vec_);
    return *this;
  }

  Vector(const Eigen::Matrix<T, Eigen::Dynamic, 1>& v):
  vec_(v)
  { }

  Vector<T>& operator=(const Vector<T>& other)
  {
    if(this != &other)
      vec_ = other.vec_; // Eigen handles deep copy

    return *this;
  }

  Vector<T>& operator*=(const T& scalar)
  {
    vec_ *= scalar;
    return *this;
  }

  Vector<T>& operator/=(const T& scalar)
  {
    if (scalar == T(0))
      throw bpp::Exception("Division by zero in Vector::operator/=");

    vec_ /= scalar;
    return *this;
  }

  bool operator==(const Vector<T>& other) const
  {
    return vec_.isApprox(other.vec_);
  }

  T operator[](size_t i) const
  {
    return eigen()(Eigen::Index(i));
  }

  T get(size_t i) const
  {
    return vec_(static_cast<Eigen::Index>(i));
  }

  void set(size_t i, const T& v)
  {
    vec_(static_cast<Eigen::Index>(i)) = v;
  }

  void setZero()
  {
    vec_.setZero();
  }

  void resize(size_t newSize)
  {
    vec_.resize(newSize);
  }

  mpfr::mpreal getMPReal(size_t index) const
  {
    return static_cast<mpfr::mpreal>(get(index));
  }

  size_t size() const
  {
    return vec_.size();
  }

  void scale(T scalar)
  {
    vec_ *= scalar;
  }

  void normalize()
  {
    T norm = vec_.norm();
    if(norm != T(0))
      vec_ /= norm;
  }

  void print(std::ostream& out) const
  {
    for(Eigen::Index i = 0; i < vec_.size(); ++i)
      out << i << ": " << vec_(i) << "\n";
  }

  std::unique_ptr<Vector<T>> clone() const
  {
    return std::make_unique<Vector<T>>(vec_);
  }

  std::unique_ptr<Vector<T>> cloneWithSize(size_t newSize) const
  {
    auto newVec = std::make_unique<Vector<T>>(newSize);
    newVec->setZero();
    return newVec;
  }

  const Eigen::Matrix<T, Eigen::Dynamic, 1>& eigen() const
  {
    return vec_;
  }

  Eigen::Matrix<T, Eigen::Dynamic, 1>& eigen()
  {
    return vec_;
  }
};

#endif
