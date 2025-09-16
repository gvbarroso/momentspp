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
#include <eigen3/Eigen/Sparse>
#include <eigen3/unsupported/Eigen/MPRealSupport>
#include <mpreal.h>
#include <vector>
#include <string>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <iostream>

template <typename T>
class Matrix
{
public:
  using Scalar = T;
  Eigen::SparseMatrix<T> mat_;

  Matrix() = default;

  Matrix(size_t rows, size_t cols) : mat_(rows, cols)
  {
  }

  Matrix(const Eigen::SparseMatrix<T>& mat) : mat_(mat)
  {
  }

  Matrix(Eigen::SparseMatrix<T>&& mat) : mat_(std::move(mat))
  {
  }

  bool operator==(const Matrix<T>& other) const
  {
    return mat_.isApprox(other.mat_);
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

  void print(const std::string& fileName) const
  {
    std::ofstream out(fileName);

    if(!out)
      throw bpp::Exception("Matrix::Failed to open file: " + fileName);

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
    result->vec_ = mat_ * vec.vec_;
    return result;
  }

  std::unique_ptr<Matrix<T>> multiply(const Matrix<T>& other) const
  {
    if(cols() != other.rows())
      throw bpp::Exception("Matrix::dimensions incompatible for multiplication!");

    Eigen::SparseMatrix<T> result = mat_ * other.mat_;
    return std::make_unique<Matrix<T>>(result);
  }

  std::unique_ptr<Matrix<T>> add(const Matrix<T>& other) const
  {
    if(rows() != other.rows() || cols() != other.cols())
      throw bpp::Exception("Matrix::dimensions must match for addition");

    Eigen::SparseMatrix<T> result = mat_ + other.mat_;
    return std::make_unique<Matrix<T>>(result);
  }

  void addInPlace(const Matrix<T>& other)
  {
    if(rows() != other.rows() || cols() != other.cols())
      throw bpp::Exception("Matrix::dimensions must match for in-place addition");

    mat_ += other.mat_;
    makeCompressed();
  }

  std::unique_ptr<Vector<T>> solve(const Vector<T>& rhs) const
  {
    if(rows() != rhs.size())
      throw bpp::Exception("Matrix and RHS dimensions do not match");

    Vector<T> result(cols());

    if constexpr(std::is_same_v<T, double>)
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
    else if constexpr(std::is_same_v<T, mpfr::mpreal>)
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

    return std::make_unique<Vector<T>>(result);
  }
};

#endif
