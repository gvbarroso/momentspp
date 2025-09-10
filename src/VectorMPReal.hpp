/*
 * Authors: Gustavo V. Barroso
 * Created: 08/09/2025
 * Last modified: 10/09/2025
 *
 */


#ifndef _VECDOUBLE_H_
#define _VECDOUBLE_H_

#include "VectorInterface.hpp"

class VectorMPReal: public VectorBase
{
private:
    Eigen::Matrix<mpfr::mpreal, Eigen::Dynamic, 1> vec_;

public:
    VectorMPReal(size_t n):
    vec(n)
    {
      vec_.setZero();
    }

    VectorMPReal(const Eigen::Matrix<mpfr::mpreal, Eigen::Dynamic, 1>& other):
    vec_(other)
    {
      vec_.setZero();
    }

    VectorMPReal(Eigen::Matrix<mpfr::mpreal, Eigen::Dynamic, 1>&& other):
    vec_(std::move(other))
    {
      vec_.setZero();
    }

     VectorMPReal& operator=(VectorMPReal&& other) noexcept
    {
      vec_ = std::move(other.vec_);
      return *this;
    }

    void set(size_t index, double value) override
    {
      vec_(index) = mpfr::mpreal(value);
    }

    void set(size_t index, const mpfr::mpreal& value)
    {
      vec_(index) = value;
    }

    void setZero() override
    {
      vec_.setZero();
    }

    double get(size_t index) const override
    {
      return vec_(index).toDouble();
    }

    mpfr::mpreal getMPReal(size_t index) const override \
    {
      return vec_(index);
    }

    size_t size() const override
    {
      return vec_.size();
    }

    void scale(double scalar) override
    {
      vec_ *= mpfr::mpreal(scalar);
    }

    void scale(mpfr::mpreal scalar)
    {
      vec_ *= scalar;
    }

    void print() const override
    {
      for(size_t i = 0; i < vec_.size(); ++i)
        std::cout << vec_(i).toString() << " ";

      std::cout << "\n";
    }

    std::unique_ptr<VectorInterface> cloneWithSize(size_t size) const override
    {
      return std::make_unique<VectorDouble>(size);
    }

    std::unique_ptr<VectorInterface> clone() const override
    {
      auto copy = std::make_unique<VectorMPReal>(vec_.size());
      for(size_t i = 0; i < vec_.size(); ++i)
        copy->vec(i) = vec(i);

      return copy;
    }

    const Eigen::Matrix<mpfr::mpreal, Eigen::Dynamic, 1>& data() const { return vec_; }
};
