/*
 * Authors: Gustavo V. Barroso
 * Created: 08/09/2025
 * Last modified: 16/09/2025
 *
 */


#ifndef _MATINTERFACE_H_
#define _MATINTERFACE_H_

#pragma once

#include <Bpp/Exceptions.h>
#include "Vector.hpp"
#include <Eigen/Sparse>
#include <unsupported/Eigen/MPRealSupport>
#include <mpreal.h>
#include <vector>
#include <string>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <iostream>

template<typename Scalar>
class Matrix
{
public:
    Eigen::SparseMatrix<Scalar> mat_;

    Matrix() = default;

    Matrix(size_t rows, size_t cols):
    mat_(rows, cols)
    { }

    Matrix(const Eigen::SparseMatrix<Scalar>& mat):
    mat_(mat)
    { }

    Matrix(Eigen::SparseMatrix<Scalar>&& mat):
    mat_(std::move(mat))
    { }

    size_t rows() const { return mat_.rows(); }
    size_t cols() const { return mat_.cols(); }

    void insert(size_t row, size_t col, Scalar val)
    {
      mat_.coeffRef(row, col) = val;
    }

    void print(const std::string& fileName) const
    {
      std::ofstream out(fileName);

      if(!out)
        throw bpp::Exception("Matrix::Failed to open file: " + fileName);

      for(int i = 0; i < mat_.outerSize(); ++i)
      {
        for(typename Eigen::SparseMatrix<Scalar>::InnerIterator it(mat_, i); it; ++it)
          out << it.row() << " " << it.col() << " " << it.value() << "\n";
      }
    }

    void zeroNegatives()
    {
      for(int i = 0; i < mat_.outerSize(); ++i)
      {
        for(typename Eigen::SparseMatrix<Scalar>::InnerIterator it(mat_, i); it; ++it)
        {
          if(it.value() < Scalar(0))
            it.valueRef() = Scalar(0);
        }
      }
    }

    void scale(Scalar scalar)
    {
      for(int k = 0; k < mat_.outerSize(); ++k)
      {
        for(typename Eigen::SparseMatrix<Scalar>::InnerIterator it(mat_, k); it; ++it)
          it.valueRef() *= scalar;
      }
    }

    void makeCompressed()
    {
      mat_.makeCompressed();
    }

    void setFromTriplets(const std::vector<Eigen::Triplet<Scalar>>& triplets)
    {
      mat_.setFromTriplets(triplets.begin(), triplets.end());
    }

    std::unique_ptr<Matrix<Scalar>> clone() const
    {
      return std::make_unique<Matrix<Scalar>>(mat_);
    }

    std::unique_ptr<Matrix<Scalar>> fetchScaled(Scalar scalar) const
    {
      auto copy = clone();
      copy->scale(scalar);
      copy->makeCompressed();
      return copy;
    }

    std::unique_ptr<Matrix<Scalar>> identity() const
    {
      std::vector<Eigen::Triplet<Scalar>> triplets;
      for(size_t i = 0; i < rows(); ++i)
        triplets.emplace_back(i, i, Scalar(1));

      auto I = std::make_unique<Matrix<Scalar>>(rows(), cols());
      I->setFromTriplets(triplets);
      I->makeCompressed();
      return I;
    }

    std::unique_ptr<Vector<Scalar>> multiply(const Vector<Scalar>& vec) const
    {
      if(cols() != vec.size())
        throw bpp::Exception("Matrix::and vector dimensions do not match");

      auto result = std::make_unique<Vector<Scalar>>(rows());
      result->vec_ = mat_ * vec.vec_;
      return result;
    }

    std::unique_ptr<Matrix<Scalar>> multiply(const Matrix<Scalar>& other) const
    {
      if(cols() != other.rows())
        throw bpp::Exception("Matrix::dimensions incompatible for multiplication");

      Eigen::SparseMatrix<Scalar> result = mat_ * other.mat_;
      return std::make_unique<Matrix<Scalar>>(result);
    }

    std::unique_ptr<Matrix<Scalar>> add(const Matrix<Scalar>& other) const
    {
      if(rows() != other.rows() || cols() != other.cols())
        throw bpp::Exception("Matrix::dimensions must match for addition");

      Eigen::SparseMatrix<Scalar> result = mat_ + other.mat_;
      return std::make_unique<Matrix<Scalar>>(result);
    }

    void addInPlace(const Matrix<Scalar>& other)
    {
      if(rows() != other.rows() || cols() != other.cols())
        throw bpp::Exception("Matrix::dimensions must match for in-place addition");

      mat_ += other.mat_;
      makeCompressed();
    }

    std::unique_ptr<Vector<Scalar>> solve(const Vector<Scalar>& rhs) const
    {
      if(rows() != rhs.size())
        throw bpp::Exception("Matrix and RHS dimensions do not match");

      Vector<Scalar> result(cols());

      if constexpr(std::is_same_v<Scalar, double>)
      {
        // Sparse LU for double
        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.compute(mat_);

        if(solver.info() != Eigen::Success)
          throw bpp::Exception("Matrix::decomposition failed");

        result.vec_ = solver.solve(rhs.vec_);

        if(solver.info() != Eigen::Success)
          throw bpp::Exception("Matrix::Solving failed");
      }

      // NOTE since matrix entries are well represented with double precision,
      // is it better (faster) to convert to double, perform Sparse LU decomposition,
      // then convert back to mpreal? May depend on condition number.
      else if constexpr(std::is_same_v<Scalar, mpfr::mpreal>)
      {
        // Dense LU for mpreal
        Eigen::Matrix<mpfr::mpreal, Eigen::Dynamic, Eigen::Dynamic> dense = mat_.toDense();
        Eigen::FullPivLU<decltype(dense)> solver;
        solver.compute(dense);

        if(!solver.isInvertible())
          throw bpp::Exception("Matrix is not invertible");

        result.vec_ = solver.solve(rhs.vec_);
      }

      else
        throw bpp::Exception("Matrix::Unsupported scalar type for solve()");

      return std::make_unique<Vector<Scalar>>(result);
    }
};


#endif
