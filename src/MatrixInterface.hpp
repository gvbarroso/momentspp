/*
 * Authors: Gustavo V. Barroso
 * Created: 08/09/2025
 * Last modified: 10/09/2025
 *
 */


#ifndef _MATINTERFACE_H_
#define _MATINTERFACE_H_

#include "eigen_mpreal_traits.hpp"

#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Core>

#include <iostream>
#include <memory>
#include <vector>
#include <mpreal.h>
#include <memory>
#include <string>
#include <fstream>

#include <Bpp/Exceptions.h>

// This is the base Matrix class for implementing runtime type polymorphism (double or mpreal)
class MatrixInterface
{
public:
    virtual ~MatrixInterface() = default;

    virtual void insert(size_t row, size_t col, double val) = 0;
    virtual void print(const std::string& fileName) const = 0;
    virtual void zeroNegatives() = 0;
    virtual void scale(double scalar) = 0;
    virtual void makeCompressed() = 0;
    virtual void setFromTriplets(const std::vector<Eigen::Triplet<double>>& triplets) = 0;
    virtual std::unique_ptr<MatrixInterface> clone() const = 0;
    virtual std::unique_ptr<MatrixInterface> identity() = 0;
    virtual std::unique_ptr<VectorInterface> multiply(const VectorInterface& vec) const = 0;
    virtual std::unique_ptr<MatrixInterface> multiply(const MatrixInterface& other) const = 0;
    virtual std::unique_ptr<MatrixInterface> add(const MatrixInterface& other) const = 0;
    virtual void addInPlace(const MatrixInterface& other) = 0;
    virtual std::unique_ptr<VectorInterface> solve(const VectorInterface& rhs) const = 0;
};

#endif
