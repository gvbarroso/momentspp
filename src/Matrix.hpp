/*
 * Authors: Gustavo V. Barroso
 * Created: 08/09/2025
 * Last modified: 25/09/2025
 *
 */

#ifndef _MATINTERFACE_H_
#define _MATINTERFACE_H_

#pragma once

#include <Bpp/Exceptions.h>
#include "Vector.hpp"
#include <eigen3/Eigen/Sparse>
#include <eigen3/unsupported/Eigen/MPRealSupport>
#include <mpreal.h>
#include <vector>
#include <string>
#include <ostream>
#include <memory>
#include <stdexcept>
#include <iostream>

template <typename T>
class Matrix
{

private:
  Eigen::SparseMatrix<T> mat_;

public:
  using Scalar = T;
  using EigenSparse = Eigen::SparseMatrix<T>;

  Matrix():
  mat_()
  { }

  Matrix(size_t rows, size_t cols):
  mat_(rows, cols)
  { }

  Matrix(const Eigen::SparseMatrix<T>& mat):
  mat_(mat)
  { }

  Matrix(Eigen::SparseMatrix<T>&& mat):
  mat_(std::move(mat))
  { }

  bool operator==(const Matrix<T>& other) const
  {
    return mat_.isApprox(other.eigen());
  }

  Matrix<T>& operator*=(const T& scalar)
  {
    mat_ *= scalar;
    return *this;
  }

  Matrix<T>& operator/=(const T& scalar)
  {
    if(scalar == T(0))
      throw bpp::Exception("Division by zero in Matrix::operator/=");

    mat_ /= scalar;
    return *this;
  }

  /// operator+= so MatrixEngine can do A += B;
  Matrix<T>& operator+=(const Matrix<T>& other)
  {
    addInPlace(other);
    return *this;
  }

  EigenSparse& eigen()
  {
    return mat_;

  }
  const EigenSparse& eigen() const
  {
    return mat_;
  }

  size_t rows() const
  {
    return mat_.rows();
  }

  size_t cols() const
  {
    return mat_.cols();
  }

  void insert(size_t row, size_t col, T val)
  {
    mat_.coeffRef(row, col) = val;
  }

 void print(std::ostream& out) const
 {
   for(int i = 0; i < mat_.outerSize(); ++i)
   {
     for(typename Eigen::SparseMatrix<T>::InnerIterator it(mat_, i); it; ++it)
       out << it.row() << " " << it.col() << " " << it.value() << "\n";
   }
 }

  void zeroNegatives()
  {
    for(int i = 0; i < mat_.outerSize(); ++i)
    {
      for(typename Eigen::SparseMatrix<T>::InnerIterator it(mat_, i); it; ++it)
      {
        if(it.value() < T(0))
          it.valueRef() = T(0);
      }
    }
  }

  void scale(T scalar)
  {
    for(int k = 0; k < mat_.outerSize(); ++k)
    {
      for(typename Eigen::SparseMatrix<T>::InnerIterator it(mat_, k); it; ++it)
        it.valueRef() *= scalar;
    }
  }

  void makeCompressed()
  {
    mat_.makeCompressed();
  }

  void prune(T threshold = T(0))
  {
    mat_.prune(threshold);
  }

  void resize(size_t rows, size_t cols)
  {
    mat_.resize(rows, cols);
  }

  void setFromTriplets(const std::vector<Eigen::Triplet<T>>& triplets)
  {
    mat_.setFromTriplets(triplets.begin(), triplets.end());
  }

  size_t nonZeros() const
  {
    return mat_.nonZeros();
  }

  bool isCompressed() const
  {
    return mat_.isCompressed();
  }

  T& coeffRef(size_t row, size_t col)
  {
    return mat_.coeffRef(row, col);
  }

  std::unique_ptr<Matrix<T>> clone() const
  {
    return std::make_unique<Matrix<T>>(mat_);
  }

  std::unique_ptr<Matrix<T>> fetchScaled(T scalar) const
  {
    auto copy = clone();
    copy->scale(scalar);
    copy->makeCompressed();
    return copy;
  }

  std::unique_ptr<Matrix<T>> identity() const
  {
    std::vector<Eigen::Triplet<T>> triplets;
    for(size_t i = 0; i < rows(); ++i)
      triplets.emplace_back(i, i, T(1));

    auto I = std::make_unique<Matrix<T>>(rows(), cols());
    I->setFromTriplets(triplets);
    I->makeCompressed();
    return I;
  }

  std::unique_ptr<Vector<T>> multiply(const Vector<T>& vec) const
  {
    if(cols() != vec.size())
       throw bpp::Exception("Matrix::Matrix and vector dimensions do not match");

    auto result = std::make_unique<Vector<T>>(rows());
    result->eigen() = mat_ * vec.eigen();
    return result;
  }

  std::unique_ptr<Vector<T>> multiply(std::unique_ptr<Vector<T>> vec) const
  {
    return multiply(*vec);
  }

  std::unique_ptr<Matrix<T>> multiply(const Matrix<T>& other) const
  {
    if(cols() != other.rows())
      throw bpp::Exception("Matrix::dimensions incompatible for multiplication!");

    Eigen::SparseMatrix<T> result = mat_ * other.eigen();
    return std::make_unique<Matrix<T>>(result);
  }

  std::unique_ptr<Matrix<T>> multiply(std::unique_ptr<Matrix<T>> other) const
  {
    return multiply(*other);
  }

  std::unique_ptr<Matrix<T>> add(const Matrix<T>& other) const
  {
    if(rows() != other.rows() || cols() != other.cols())
      throw bpp::Exception("Matrix::dimensions must match for addition");

    Eigen::SparseMatrix<T> result = mat_ + other.eigen();
    return std::make_unique<Matrix<T>>(result);
  }

  void addInPlace(const Matrix<T>& other)
  {
    if(rows() != other.rows() || cols() != other.cols())
      throw bpp::Exception("Matrix::dimensions must match for in-place addition");

    mat_ += other.eigen();
    makeCompressed();
  }

  std::unique_ptr<Vector<T>> solve(const Vector<T>& rhs) const
  {
    if(cols() != rhs.size())
      throw bpp::Exception("Matrix::Matrix and RHS dimensions do not match");

    Vector<T> result(rows());

    if constexpr(std::is_same_v<T, double>)
    {
      Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
      const Eigen::SparseMatrix<double>& A = eigen();
      solver.compute(A);

      if(solver.info() != Eigen::Success)
        throw bpp::Exception("Matrix::decomposition failed");

      result.eigen() = solver.solve(rhs.eigen());

      if(solver.info() != Eigen::Success)
        throw bpp::Exception("Matrix::Solving failed");
    }

    else if constexpr(std::is_same_v<T, mpfr::mpreal>)
    {
      // convert matrix to double to use SparseLU
      const Eigen::SparseMatrix<mpfr::mpreal>& A_mp = eigen();
      Eigen::SparseMatrix<double> A_d(A_mp.rows(), A_mp.cols());
      A_d.reserve(A_mp.nonZeros());

      for(int k = 0; k < A_mp.outerSize(); ++k)
      {
        for(typename Eigen::SparseMatrix<mpfr::mpreal>::InnerIterator it(A_mp, k); it; ++it)
          A_d.coeffRef(it.row(), it.col()) = static_cast<double>(it.value());
      }

      A_d.makeCompressed();

      // convert RHS to double
      const Eigen::Matrix<mpfr::mpreal, Eigen::Dynamic, 1>& rhs_mp = rhs.eigen();
      Eigen::VectorXd rhs_d(rhs_mp.size());

      for(Eigen::Index i = 0; i < rhs_mp.size(); ++i)
        rhs_d(i) = static_cast<double>(rhs_mp(i));

      // solve in double
      Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
      solver.compute(A_d);

      if(solver.info() != Eigen::Success)
        throw bpp::Exception("Matrix::decomposition failed");

      Eigen::VectorXd x_d = solver.solve(rhs_d);

      if(solver.info() != Eigen::Success)
        throw bpp::Exception("Matrix::Solving failed");

      // convert solution back to mpfr
      for(Eigen::Index i = 0; i < x_d.size(); ++i)
        result.eigen()(i) = mpfr::mpreal(x_d(i));
    }

    else
      throw bpp::Exception("Matrix::Unsupported scalar type for solve()");

    return std::make_unique<Vector<T>>(result);
  }

  // Add identity: M ← I + M
  void addIdentity()
  {
    for(size_t i = 0; i < rows(); ++i)
      mat_.coeffRef(i, i) += T(1);
    makeCompressed();
  }

  // Alias for makeCompressed (named compress in MatrixEngine)
  void compress()
  {
    makeCompressed();
  }

};

#endif
