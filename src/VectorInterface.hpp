/*
 * Authors: Gustavo V. Barroso
 * Created: 08/09/2025
 * Last modified: 10/09/2025
 *
 */

#ifndef _VECINTERFACE_H_
#define _VECINTERFACE_H_

#include "eigen_mpreal_traits.hpp"

#include <Eigen/Sparse>
#include <iostream>
#include <memory>
#include <vector>
#include <mpreal.h>
#include <memory>


class VectorInterface
{
public:
    virtual ~VectorInterface() = default;

    virtual void set(size_t index, double value) = 0;
    virtual void setZero() = 0;
    virtual double get(size_t index) const = 0;
    virtual mpfr::mpreal getMPReal(size_t index) const { return get(index); } // returns double by default
    virtual size_t size() const = 0;
    virtual void scale(double scalar) = 0;
    virtual void print() const = 0;
    virtual std::unique_ptr<VectorInterface> clone() const = 0;
    virtual std::unique_ptr<VectorInterface> cloneWithSize(size_t size) const = 0;
};

#endif
