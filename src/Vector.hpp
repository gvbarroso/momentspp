/*
 * Authors: Gustavo V. Barroso
 * Created: 08/09/2025
 * Last modified: 15/09/2025
 *
 */

#ifndef _VECINTERFACE_H_
#define _VECINTERFACE_H_

#pragma once

#include <Eigen/Core>
#include <unsupported/Eigen/MPRealSupport>
#include <mpreal.h>
#include <memory>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdexcept>

template<typename Scalar>
class Vector
{
public:
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> vec_;

    Vector() = default;

    explicit Vector(size_t size):
    vec_(size)
    {
      setZero();
    }

    Vector(const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& v):
    vec_(v)
    {
      setZero();
    }

    Vector(Eigen::Matrix<Scalar, Eigen::Dynamic, 1>&& other):
    vec_(std::move(other))
    {
      setZero();
    }

    Vector<Scalar>& operator=(const Vector<Scalar>& other)
    {
      if(this != &other)
        vec_ = other.vec_;  // Eigen handles deep copy

      return *this;
    }

    void set(size_t index, Scalar value)
    {
      if(index >= vec_.size())
        throw std::out_of_range("Vector index out of bounds");

      vec_(index) = value;
    }

    void setZero()
    {
      vec_.setZero();
    }

    Scalar get(size_t index) const
    {
      if(index >= vec_.size())
        throw std::out_of_range("Vector index out of bounds");

      return vec_(index);
    }

    mpfr::mpreal getMPReal(size_t index) const
    {
      return static_cast<mpfr::mpreal>(get(index));
    }

    size_t size() const
    {
      return vec_.size();
    }

    void scale(Scalar scalar)
    {
      vec_ *= scalar;
    }

    void print() const
    {
      for(size_t i = 0; i < vec_.size(); ++i)
        std::cout << i << ": " << vec_(i) << "\n";
    }

    std::unique_ptr<Vector<Scalar>> clone() const
    {
      return std::make_unique<Vector<Scalar>>(vec_);
    }

    std::unique_ptr<Vector<Scalar>> cloneWithSize(size_t newSize) const
    {
      auto newVec = std::make_unique<Vector<Scalar>>(newSize);
      newVec->setZero();
      return newVec;
    }

    const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& eigen() const { return vec_; }
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& eigen() { return vec_; }
};


#endif
