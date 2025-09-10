/*
 * Authors: Gustavo V. Barroso
 * Created: 08/09/2025
 * Last modified:  10/09/2025
 */


#ifndef _VECDOUBLE_H_
#define _VECDOUBLE_H_

#include "VectorInterface.hpp"

class VectorDouble : public VectorInterface
{
private:
    Eigen::VectorXd vec_;

public:
    VectorDouble(size_t n):
    vec_(n)
    {
      vec_.setZero();
    }

    VectorDouble(const Eigen::Matrix<double, Eigen::Dynamic, 1>& other):
    vec_(other)
    {
      vec_.setZero();
    }

    VectorDouble(Eigen::Matrix<double, Eigen::Dynamic, 1>&& other):
    vec_(std::move(other))
    {
      vec_.setZero();
    }

    /*
    VectorDouble& operator=(VectorDouble&& other) noexcept
    {
      vec_ = std::move(other.vec_);
      return *this;
    }
    */

    void set(size_t index, double value) override
    {
      vec_(index) = value;
    }

    void setZero() override
    {
      vec_.setZero();
    }

    double get(size_t index) const override
    {
      return vec_(index);
    }

    size_t size() const override
    {
      return vec_.size();
    }

    void scale(double scalar) override {
      vec_ *= scalar;
    }

    void print() const override
    {
      std::cout << vec_.transpose() << "\n";
    }

    std::unique_ptr<VectorInterface> cloneWithSize(size_t size) const override
    {
      return std::make_unique<VectorDouble>(size);
    }

    std::unique_ptr<VectorInterface> clone() const override
    {
      auto copy = std::make_unique<VectorDouble>(vec_.size());
      for(size_t i = 0; i < vec_.size(); ++i)
        copy->set(i, vec_(i));

      return copy;
    }

    const Eigen::VectorXd& data() const { return vec_; }
};

#endif
