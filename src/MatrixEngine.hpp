/*
 * Authors: Gustavo V. Barroso
 * Created: 12/09/2025
 * Last modified: 16/09/2025
 *
 */

#pragma once

#include <mpreal.h>

#include "Matrix.hpp"
#include "Vector.hpp"

#include <variant>
#include <string>
#include <memory>
#include <stdexcept>
#include <vector>
#include <eigen3/Eigen/Sparse>

class MatrixEngine
{
public:
  using MatrixVariant = std::variant<Matrix<double>, Matrix<mpfr::mpreal>>;
  using VectorVariant = std::variant<Vector<double>, Vector<mpfr::mpreal>>;

  MatrixVariant matrix;
  VectorVariant vector;
  bool useMPRealFlag;

  MatrixEngine(bool useMPRealPrecision = false):
  useMPRealFlag(useMPRealPrecision)
  { }

  bool useMPReal() const
  {
    return useMPRealFlag;
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
    std::visit([](auto& vec)
               { vec.setZero(); },
               vector);
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
        for(int j = 0; j < mat.cols(); ++j)
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
    std::visit([](auto& vec)
    { vec.normalize();
    }, vector);
  }

  void logMatrixStats(const std::string& label)
  {
    std::visit([&](const auto& mat)
    {
      std::cout << label << ": " << mat.rows() << "x" << mat.cols() << ", nonzeros = " << mat.nonZeros() << "\n";
    }, matrix);
  }

  using VectorVariantEigen = std::variant<Eigen::VectorXd, Eigen::Matrix<mpfr::mpreal, Eigen::Dynamic, 1>>;
  VectorVariantEigen getRawVector() const
  {
    return std::visit([](const auto& vec) -> VectorVariantEigen
    { return vec.vec_;
    }, vector);
  }

  using SparseMatrixVariant = std::variant<Eigen::SparseMatrix<double>, Eigen::SparseMatrix<mpfr::mpreal>>;
  SparseMatrixVariant getRawMatrix() const
  {
    return std::visit([](const auto& mat) -> SparseMatrixVariant
        { return mat.mat_; },
                      matrix);
  }

  template <typename Scalar>
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
    if(useMPReal())
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
    std::visit( [&](auto& mat)
    {
      using Scalar = typename std::decay_t<decltype(mat)>::Scalar;
      std::vector<Eigen::Triplet<Scalar>> converted;

      for(const auto& t : triplets)
        converted.emplace_back(t.row(), t.col(), static_cast<Scalar>(t.value()));

      mat.setFromTriplets(converted);
    }, matrix);
  }

  size_t matrixRows() const
  {
    return std::visit([](const auto& mat)
    { return mat.rows();
    }, matrix);
  }

  size_t matrixCols() const
  {
    return std::visit([](const auto& mat)
    { return mat.cols();
    }, matrix);
  }

  size_t vectorSize() const
  {
    return std::visit([](const auto& vec)
    { return vec.size();
    }, vector);
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
    { mat.makeCompressed();
    }, matrix);
  }

  void compressInPlace()
  {
    std::visit([](auto& mat)
    { mat.makeCompressed();
    }, matrix);
  }

  void pruneInPlace(double threshold = 0.)
  {
    std::visit([threshold](auto& mat)
    { mat.prune(static_cast<typename std::decay_t<decltype(mat)>::Scalar>(threshold));
    }, matrix);
  }

  void zeroMatrixNegatives()
  {
    std::visit([](auto& mat)
    { mat.zeroNegatives();
    }, matrix);
  }

  void printMatrix(const std::string& filename)
  {
    std::visit([&](auto& mat)
    { mat.print(filename);
    }, matrix);
  }

  void addMatrixInPlace(const MatrixVariant& other)
  {
    if(matrix.index() != other.index())
      throw bpp::Exception("MatrixEngine::types do not match for in-place addition.");

    if(std::holds_alternative<Matrix<double>>(matrix))
    {
      auto& mat = std::get<Matrix<double>>(matrix);
      const auto& otherMat = std::get<Matrix<double>>(other);
      mat.addInPlace(otherMat);
    }

    else
    {
      auto& mat = std::get<Matrix<mpfr::mpreal>>(matrix);
      const auto& otherMat = std::get<Matrix<mpfr::mpreal>>(other);
      mat.addInPlace(otherMat);
    }
  }

  void addIdentityInPlace()
  {
    std::visit([](auto& mat)
    {
      using Scalar = typename std::decay_t<decltype(mat)>::Scalar;
      const Eigen::Index size = mat.rows();

      for(Eigen::Index i = 0; i < size; ++i)
        mat.mat_.coeffRef(i, i) += Scalar(1);
    }, matrix);
  }

  MatrixVariant addMatrix(const MatrixVariant& other)
  {
    if(matrix.index() != other.index())
      throw bpp::Exception("MatrixEngine::types do not match for addition.");

    if(std::holds_alternative<Matrix<double>>(matrix))
    {
      const auto& mat = std::get<Matrix<double>>(matrix);
      const auto& otherMat = std::get<Matrix<double>>(other);
      return *mat.add(otherMat);
    }

    else
    {
      const auto& mat = std::get<Matrix<mpfr::mpreal>>(matrix);
      const auto& otherMat = std::get<Matrix<mpfr::mpreal>>(other);
      return *mat.add(otherMat);
    }
  }

  MatrixVariant multiplyMatrix(const MatrixVariant& other)
  {
    if(matrix.index() != other.index())
      throw bpp::Exception("MatrixEngine::types do not match for multiplication.");

    if(std::holds_alternative<Matrix<double>>(matrix))
    {
      const auto& mat = std::get<Matrix<double>>(matrix);
      const auto& otherMat = std::get<Matrix<double>>(other);
      return *mat.multiply(otherMat);
    }

    else
    {
      const auto& mat = std::get<Matrix<mpfr::mpreal>>(matrix);
      const auto& otherMat = std::get<Matrix<mpfr::mpreal>>(other);
      return *mat.multiply(otherMat);
    }
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
    { return *mat.identity();
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
    { vec.setZero();
    }, vector);
  }

  void scaleVector(double scalar)
  {
    std::visit( [=](auto& vec)
    {
      using Scalar = typename std::decay_t<decltype(vec)>::Scalar;
      vec.scale(static_cast<Scalar>(scalar));
    }, vector);
  }

  void printVector()
  {
    std::visit([](auto& vec)
    { vec.print();
    }, vector);
  }

  VectorVariant cloneVector()
  {
    return std::visit([](const auto& vec) -> VectorVariant
    { return *vec.clone();
    }, vector);
  }

  VectorVariant cloneVectorWithSize(size_t size)
  {
    return std::visit([=](const auto& vec) -> VectorVariant
    { return *vec.cloneWithSize(size);
    }, vector);
  }

  // Matrix-vector operations
  VectorVariant multiplyMatrixVector() const
  {
    return std::visit([](const auto& mat, const auto& vec) -> VectorVariant
    {
      using MatType = std::decay_t<decltype(mat)>;
      using VecType = std::decay_t<decltype(vec)>;

      if constexpr (std::is_same<typename MatType::Scalar, typename VecType::Scalar>::value)
        return VectorVariant{*mat.multiply(vec)};

      else
        throw bpp::Exception("MatrixEngine::multiplyMatrixVector: mismatched scalar types between matrix and vector.");
    }, matrix, vector);
  }

  VectorVariant solveSystem() const
  {
    return std::visit([](const auto& mat, const auto& vec) -> VectorVariant
    {
      using MatType = std::decay_t<decltype(mat)>;
      using VecType = std::decay_t<decltype(vec)>;

      if constexpr (std::is_same<typename MatType::Scalar, typename VecType::Scalar>::value)
        return VectorVariant{*mat.solve(vec)};

      else
        throw bpp::Exception("MatrixEngine::solveSystem: mismatched scalar types between matrix and vector.");
    }, matrix, vector);
  }
};
