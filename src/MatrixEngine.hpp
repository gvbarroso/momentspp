/*
 * Authors: Gustavo V. Barroso
 * Created: 12/09/2025
 * Last modified: 15/09/2025
 *
 */

#pragma once

#include "Matrix.hpp"
#include "Vector.hpp"

#include <variant>
#include <string>
#include <memory>
#include <stdexcept>
#include <vector>
#include <Eigen/Sparse>

class MatrixEngine
{
public:
  using MatrixVariant = std::variant<Matrix<double>, Matrix<mpfr::mpreal>>;
  using VectorVariant = std::variant<Vector<double>, Vector<mpfr::mpreal>>;

  MatrixVariant matrix;
  VectorVariant vector;
  bool useMPReal;

  MatrixEngine(bool useMPRealPrecision = false):
  useMPReal(useMPRealPrecision)
  { }

  bool useMPReal() const
  {
    return useMPReal;
  }

  void setMatrix(const MatrixVariant& newMatrix)
  {
    matrix = newMatrix;
  }

  void setVector(const VectorVariant& newVector)
  {
    vector = newVector;
  }

  void resetVector()
  {
    std::visit([](auto& vec) {
      vec.setZero();
    }, vector);
  }

  void setMatrixFromEigen(const Eigen::SparseMatrix<double>& mat)
  {
    std::visit([&](auto& m)
    {
      using Scalar = typename std::decay_t<decltype(m)>::Scalar;
      Eigen::SparseMatrix<Scalar> converted(mat.rows(), mat.cols());

      std::vector<Eigen::Triplet<Scalar>> triplets;
      for(int k = 0; k < mat.outerSize(); ++k)
      {
        for(Eigen::SparseMatrix<double>::InnerIterator it(mat, k); it; ++it)
          triplets.emplace_back(it.row(), it.col(), static_cast<Scalar>(it.value()));
      }

      converted.setFromTriplets(triplets.begin(), triplets.end());
      m.mat_ = converted;
    }, matrix);
  }

  void setMatrixFromDense(const Eigen::MatrixXd& mat)
  {
    std::visit([&](auto& m)
    {
      using Scalar = typename std::decay_t<decltype(m)>::Scalar;
      Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> converted(mat.rows(), mat.cols());

      for(int i = 0; i < mat.rows(); ++i)
      {
        for (int j = 0; j < mat.cols(); ++j)
          converted(i, j) = static_cast<Scalar>(mat(i, j));
      }

      m.mat_ = converted.sparseView(); // convert to sparse if needed
    }, matrix);
  }

  void setVectorFromEigen(const Eigen::VectorXd& vec)
  {
    std::visit([&](auto& v)
    {
      using Scalar = typename std::decay_t<decltype(v)>::Scalar;
      Eigen::Matrix<Scalar, Eigen::Dynamic, 1> converted(vec.size());

      for(int i = 0; i < vec.size(); ++i)
        converted(i) = static_cast<Scalar>(vec(i));

      v.vec_ = converted;
    }, vector);
  }

  void setVectorFromMPReal(const Eigen::Matrix<mpfr::mpreal, Eigen::Dynamic, 1>& vec)
  {
    std::visit([&](auto& v)
    {
      using Scalar = typename std::decay_t<decltype(v)>::Scalar;
      Eigen::Matrix<Scalar, Eigen::Dynamic, 1> converted(vec.size());

      for(int i = 0; i < vec.size(); ++i)
        converted(i) = static_cast<Scalar>(vec(i));

      v.vec_ = converted;
    }, vector);
  }

  MatrixVariant& getMatrixVariant()
  {
    return matrix;
  }

  const MatrixVariant& getMatrixVariant() const
  {
    return matrix;
  }

  VectorVariant& getVectorVariant()
  {
    return vector;
  }

  const VectorVariant& getVectorVariant() const
  {
    return vector;
  }

  void normalizeVector()
  {
    std::visit([](auto& vec) {
      vec.normalize();
    }, vector);
  }

  void logMatrixStats(const std::string& label)
  {
    std::visit([&](const auto& mat) {
      std::cout << label << ": " << mat.rows() << "x" << mat.cols()
                << ", nonzeros = " << mat.nonZeros() << "\n";
    }, matrix);
  }

  using VectorVariantEigen = std::variant< Eigen::VectorXd, Eigen::Matrix<mpfr::mpreal, Eigen::Dynamic, 1>>;
  VectorVariantEigen getRawVector() const
  {
    return std::visit([](const auto& vec)
    {
      return vec.vec_;
    }, vector);
  }

  using SparseMatrixVariant = std::variant< Eigen::SparseMatrix<double>, Eigen::SparseMatrix<mpfr::mpreal>>;
  SparseMatrixVariant getRawMatrix() const
  {
    return std::visit([](const auto& mat)
    {
      return mat.mat_;
    }, matrix);
  }

  template<typename Scalar>
  Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> getMatrixAsDense() const
  {
    return std::visit([](const auto& mat)
    {
      using MatScalar = typename std::decay_t<decltype(mat)>::Scalar;
      static_assert(std::is_same_v<MatScalar, Scalar>, "Requested type does not match stored matrix type.");
      return Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>(mat.mat_);
    }, matrix);
  }

  void initialize(size_t rows, size_t cols, size_t vecSize)
  {
    if(useMPReal)
    {
      matrix = Matrix<mpfr::mpreal>(rows, cols);
      vector = Vector<mpfr::mpreal>(vecSize);
    }

    else
    {
      matrix = Matrix<double>(rows, cols);
      vector = Vector<double>(vecSize);
    }
  }

  void insertMatrixValue(size_t row, size_t col, double val)
  {
    std::visit([=](auto& mat)
    {
      using Scalar = typename std::decay_t<decltype(mat)>::Scalar;
      mat.insert(row, col, static_cast<Scalar>(val));
    }, matrix);
  }

  void setMatrixFromTriplets(const std::vector<Eigen::Triplet<double>>& triplets)
  {
    std::visit([&](auto& mat)
    {
      using Scalar = typename std::decay_t<decltype(mat)>::Scalar;
      std::vector<Eigen::Triplet<Scalar>> converted;
      for (const auto& t : triplets)
      {
        converted.emplace_back(t.row(), t.col(), static_cast<Scalar>(t.value()));
      }
      mat.setFromTriplets(converted);
    }, matrix);
  }

  void scaleMatrix(double scalar)
  {
    std::visit([=](auto& mat)
    {
      using Scalar = typename std::decay_t<decltype(mat)>::Scalar;
      mat.scale(static_cast<Scalar>(scalar));
    }, matrix);
  }

  void makeMatrixCompressed()
  {
    std::visit([](auto& mat)
    {
      mat.makeCompressed();
    }, matrix);
  }

  void compressInPlace()
  {
    std::visit([](auto& mat) {
      mat.makeCompressed();
    }, matrix_);
  }

  void pruneInPlace(double threshold = 0.)
  {
    std::visit([threshold](auto& mat) {
      mat.prune(static_cast<typename std::decay_t<decltype(mat)>::Scalar>(threshold));
    }, matrix_);
  }

  void zeroMatrixNegatives()
  {
    std::visit([](auto& mat)
    {
      mat.zeroNegatives();
    }, matrix);
  }

  void printMatrix(const std::string& filename)
  {
    std::visit([&](auto& mat)
    {
      mat.print(filename);
    }, matrix);
  }

  void addMatrixInPlace(const MatrixVariant& other)
  {
    std::visit([&](auto& mat, const auto& otherMat)
    {
      mat.addInPlace(otherMat);
    }, matrix, other);
  }

  void addIdentityInPlace()
  {
    std::visit([](auto& mat) {
      using Scalar = typename std::decay_t<decltype(mat)>::Scalar;
      const Eigen::Index size = mat.rows();

      for(Eigen::Index i = 0; i < size; ++i)
        mat.coeffRef(i, i) += Scalar(1);
    }, matrix_);
  }

  MatrixVariant addMatrix(const MatrixVariant& other)
  {
    return std::visit([](const auto& mat, const auto& otherMat) -> MatrixVariant
    {
      return *mat.add(otherMat);
    }, matrix, other);
  }

  MatrixVariant multiplyMatrix(const MatrixVariant& other)
  {
    return std::visit([](const auto& mat, const auto& otherMat) -> MatrixVariant
    {
      return *mat.multiply(otherMat);
    }, matrix, other);
  }

  MatrixVariant fetchScaledMatrix(double scalar)
  {
    return std::visit([=](const auto& mat) -> MatrixVariant
    {
      using Scalar = typename std::decay_t<decltype(mat)>::Scalar;
      return *mat.fetchScaled(static_cast<Scalar>(scalar));
    }, matrix);
  }

  MatrixVariant identityMatrix()
  {
    return std::visit([](const auto& mat) -> MatrixVariant
    {
      return *mat.identity();
    }, matrix);
  }

  void setVectorValue(size_t index, double val)
  {
    std::visit([=](auto& vec)
    {
      using Scalar = typename std::decay_t<decltype(vec)>::Scalar;
      vec.set(index, static_cast<Scalar>(val));
    }, vector);
  }

  void setVectorZero()
  {
    std::visit([](auto& vec)
    {
      vec.setZero();
    }, vector);
  }

  void scaleVector(double scalar)
  {
    std::visit([=](auto& vec)
    {
      using Scalar = typename std::decay_t<decltype(vec)>::Scalar;
      vec.scale(static_cast<Scalar>(scalar));
    }, vector);
  }

  void printVector()
  {
    std::visit([](auto& vec)
    {
      vec.print();
    }, vector);
  }

  VectorVariant cloneVector()
  {
    return std::visit([](const auto& vec) -> VectorVariant
    {
      return *vec.clone();
    }, vector);
  }

  VectorVariant cloneVectorWithSize(size_t size)
  {
    return std::visit([=](const auto& vec) -> VectorVariant
    {
      return *vec.cloneWithSize(size);
    }, vector);
  }

  // Matrix-vector operations
  VectorVariant multiplyMatrixVector()
  {
    return std::visit([](const auto& mat, const auto& vec) -> VectorVariant
    {
      return *mat.multiply(vec);
    }, matrix, vector);
  }

  VectorVariant solveSystem()
  {
    return std::visit([](const auto& mat, const auto& vec) -> VectorVariant
    {
      return *mat.solve(vec);
    }, matrix, vector);
  }
};
