/*
 * Authors: Gustavo V. Barroso
 * Created: 12/09/2025
 * Last modified: 18/09/2025
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

  using MatrixVariantEigen = std::variant<Eigen::SparseMatrix<double>, Eigen::SparseMatrix<mpfr::mpreal>>;
  using VectorVariantEigen = std::variant<Eigen::VectorXd, Eigen::Matrix<mpfr::mpreal, Eigen::Dynamic, 1>>;

  MatrixVariant matrix;
  VectorVariant vector;
  bool useMPRealFlag;

  MatrixEngine(bool useMPRealPrecision = false):
  useMPRealFlag(useMPRealPrecision)
  { }

  MatrixEngine(std::variant<Matrix<double>, Matrix<mpfr::mpreal>> mat)
  {
    matrix = std::move(mat);
  }

  MatrixEngine& operator*=(double scalar)
  {
    std::visit([scalar](auto& mat)
    {
      mat *= scalar;
    }, matrix);

    std::visit([scalar](auto& vec)
    {
      vec *= scalar;
    }, vector);

    return *this;
  }

  MatrixEngine& operator/=(double scalar)
  {
    if(scalar == 0.0)
      throw bpp::Exception("Division by zero in MatrixEngine::operator/=");

    std::visit([scalar](auto& mat)
    {
      mat /= scalar;
    }, matrix);

    std::visit([scalar](auto& vec)
    {
      vec /= scalar;
    }, vector);

    return *this;
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

  bool useMPReal() const
  {
    return useMPRealFlag;
  }

  VectorVariantEigen toEigenVectorVariant() const
  {
    return std::visit([](const auto& vec) -> VectorVariantEigen
    {
      return VectorVariantEigen(vec.eigen());
    }, vector);
  }

  MatrixVariantEigen toEigenMatrixVariant() const
  {
    return std::visit([](const auto& mat) -> MatrixVariantEigen {
      return MatrixVariantEigen(mat.eigen());
    }, matrix);
  }

  void setMatrix(const MatrixVariant& newMatrix)
  {
    matrix = newMatrix;
  }

  void setMatrix(std::unique_ptr<MatrixVariant> newMatrix)
  {
    setMatrix(*newMatrix);
  }

  void setVector(const VectorVariant& newVector)
  {
    vector = newVector;
  }

  void setVector(std::unique_ptr<VectorVariant> newVector)
  {
    vector = std::move(*newVector);
  }

  void resetVector()
  {
    std::visit([](auto& vec)
               { vec.setZero(); },
               vector);
  }

  void setMatrixFromEigen(const MatrixVariantEigen& eigenVar)
  {
    std::visit([&](auto const& matEigen)
    {
      using EigenMatT = std::decay_t<decltype(matEigen)>;
      using Scalar    = typename EigenMatT::Scalar;

      // Construct your Matrix<Scalar> wrapper directly from the Eigen sparse matrix
      Matrix<Scalar> wrapper(matEigen);

      // Properly construct the internal variant by specifying in_place_type
      MatrixVariant mv(std::in_place_type<Matrix<Scalar>>, std::move(wrapper));

      // Call the existing setter
      setMatrix(mv);
    }, eigenVar);
  }

  void setMatrixFromEigen(MatrixVariantEigen&& eigenVar)
  {
    // Forward to the const-lvalue overload
    setMatrixFromEigen(static_cast<const MatrixVariantEigen&>(eigenVar));
  }


  void setMatrixFromDense(const Eigen::MatrixXd& mat)
  {
    std::visit([&](auto& m) -> void
    {
      using WrapT = std::decay_t<decltype(m)>;
      using Scalar = typename WrapT::Scalar;

      Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> denseConv(mat.rows(), mat.cols());
      for(int i = 0; i < mat.rows(); ++i)
      {
        for(int j = 0; j < mat.cols(); ++j)
          denseConv(i, j) = static_cast<Scalar>(mat(i, j));
      }

      m.eigen() = denseConv.sparseView();
    }, matrix);
  }

  void setVectorFromEigen(const VectorVariantEigen& eigenVar)
  {
    std::visit([&](auto const& vecEigen)
    {
      using EigenVecT = std::decay_t<decltype(vecEigen)>;
      using Scalar    = typename EigenVecT::Scalar;

      // Build Vector<Scalar> wrapper and assign the Eigen vector
      Vector<Scalar> wrapper(static_cast<size_t>(vecEigen.size()));
      wrapper.eigen() = vecEigen;

      // Construct the internal variant correctly
      VectorVariant vv(std::in_place_type<Vector<Scalar>>, std::move(wrapper));
      setVector(vv);
    }, eigenVar);
  }

  void setVectorFromEigen(VectorVariantEigen&& eigenVar)
  {
    setVectorFromEigen(static_cast<const VectorVariantEigen&>(eigenVar));
  }

  void setVectorFromMPReal(const Eigen::Matrix<mpfr::mpreal, Eigen::Dynamic, 1>& vec)
  {
    std::visit([&](auto& v) -> void
    {
      using WrapT = std::decay_t<decltype(v)>;
      using Scalar = typename WrapT::Scalar;

      Eigen::Matrix<Scalar, Eigen::Dynamic, 1> converted(vec.size());
      for(int i = 0; i < vec.size(); ++i)
        converted(i) = static_cast<Scalar>(vec(i));

      v.eigen() = std::move(converted);
    }, vector);
  }

  std::unique_ptr<MatrixEngine> clone() const
  {
    auto copy = std::make_unique<MatrixEngine>(useMPRealFlag);

    copy->matrix = std::visit([](const auto& mat) -> MatrixVariant
    {
        return MatrixVariant{*mat.clone()};
    }, matrix);

    copy->vector = std::visit([](const auto& vec) -> VectorVariant
    {
        return VectorVariant{*vec.clone()};
    }, vector);

    return copy;
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

  template <typename Scalar>
  Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> getMatrixAsDense() const
  {
    return std::visit([](const auto& mat)
    {
      using MatScalar = typename std::decay_t<decltype(mat)>::Scalar;
      static_assert(std::is_same_v<MatScalar, Scalar>, "Requested type does not match stored matrix type.");
      return Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>(mat.eigen());
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
        mat.eigen().coeffRef(i, i) += Scalar(1);
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
