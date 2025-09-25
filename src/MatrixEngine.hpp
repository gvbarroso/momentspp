/*
 * Authors: Gustavo V. Barroso
 * Created: 08/09/2025
 * Last modified: 24/09/2025
 *
 */


#ifndef _MATRIXENGINE_HPP_
#define _MATRIXENGINE_HPP_

#include "VariantUtils.hpp"
#include "Matrix.hpp" //  Matrix<Scalar> wrapper
#include "Vector.hpp" // Vector<Scalar> wrapper

#include <variant>
#include <memory>
#include <Bpp/Exceptions.h>

class MatrixEngine
{
public:
  //----------------------------------------------------------------------
  // Wrapper types (thin wrappers around Eigen)
  //----------------------------------------------------------------------
  using DoubleMatrixWrap = Matrix<double>;
  using MPRealMatrixWrap = Matrix<mpfr::mpreal>;
  using DoubleVectorWrap = Vector<double>;
  using MPRealVectorWrap = Vector<mpfr::mpreal>;

  //----------------------------------------------------------------------
  // Pure‐Eigen types (used by Epoch for integrate, print, etc.)
  //----------------------------------------------------------------------
  using DoubleMatrixEigen = Eigen::SparseMatrix<double>;
  using MPRealMatrixEigen = Eigen::SparseMatrix<mpfr::mpreal>;
  using DoubleVectorEigen = Eigen::Matrix<double, Eigen::Dynamic, 1>;
  using MPRealVectorEigen = Eigen::Matrix<mpfr::mpreal, Eigen::Dynamic, 1>;

  //----------------------------------------------------------------------
  // Variants
  //----------------------------------------------------------------------
  using MatrixVariant = std::variant<DoubleMatrixWrap, MPRealMatrixWrap>;
  using VectorVariant = std::variant<DoubleVectorWrap, MPRealVectorWrap>;
  using MatrixEigenVariant = std::variant<DoubleMatrixEigen, MPRealMatrixEigen>;
  using VectorEigenVariant = std::variant<DoubleVectorEigen, MPRealVectorEigen>;

private:
  MatrixVariant matWrap_;
  VectorVariant vecWrap_;

  int rows_{0};
  int cols_{0};

public:
  explicit MatrixEngine(bool useMPReal = false)
  {
    if(useMPReal)
      matWrap_ = MPRealMatrixWrap{};

    else
      matWrap_ = DoubleMatrixWrap{};

    std::visit(overloaded{ [&](auto const &M)
    {
      rows_ = M.rows();
      cols_ = M.cols();
    }
    }, matWrap_);
  }

  MatrixEngine(const MatrixVariant& Mwrap, const VectorVariant& Vwrap = VectorVariant{}):
  matWrap_(Mwrap),
  vecWrap_(Vwrap)
  {
    std::visit(overloaded{[&](auto const &M)
    {
      rows_ = M.rows();
      cols_ = M.cols();
    }
    }, matWrap_);
  }

  MatrixEngine(const MatrixEngine&) = default;
  MatrixEngine(MatrixEngine&&) noexcept = default;
  MatrixEngine& operator=(const MatrixEngine&) = default;
  MatrixEngine& operator=(MatrixEngine&&) noexcept = default;

  ~MatrixEngine() = default;

  std::unique_ptr<MatrixEngine> clone() const
  {
    return std::make_unique<MatrixEngine>(*this);
  }

  MatrixVariant&& moveMatrixVariant()
  {
    return std::move(matWrap_);
  }

  VectorVariant&& moveVectorVariant()
  {
    return std::move(vecWrap_);
  }

  //----------------------------------------------------------------------
  // Accessors for wrapper variants (used by AbstractOperator)
  //----------------------------------------------------------------------
  const MatrixVariant& getMatrixVariant() const
  {
    return matWrap_;
  }

  const VectorVariant& getVectorVariant() const
  {
    return vecWrap_;
  }

  //----------------------------------------------------------------------
  // Expose pure‐Eigen variants (used by Epoch::integrate, printTransitionMat)
  //----------------------------------------------------------------------
  MatrixEigenVariant toEigenMatrixVariant() const
  {
    MatrixEigenVariant out;
    std::visit(overloaded
    {
      [&](DoubleMatrixWrap const &W) { out = W.eigen(); },
      [&](MPRealMatrixWrap const &W) { out = W.eigen(); }
    }, matWrap_);
    return out;
  }

  VectorEigenVariant toEigenVectorVariant() const
  {
    VectorEigenVariant out;
    std::visit(overloaded
    {
      [&](DoubleVectorWrap const &W) { out = W.eigen(); },
      [&](MPRealVectorWrap const &W) { out = W.eigen(); }
    }, vecWrap_);
    return out;
  }

  void setVector(VectorVariant&& v)
  {
    vecWrap_ = std::move(v);
  }

  void setMatrix(MatrixVariant&& m)
  {
    matWrap_ = std::move(m);
  }

  //std::unique_ptr<MatrixEngine> me = MatrixEngine::createEmpty<double>(size);
  //std::unique_ptr<MatrixEngine> me = MatrixEngine::createEmpty<mpfr::mpreal>(size)
  template <typename Scalar>
  static std::unique_ptr<MatrixEngine> createEmpty(size_t dim) // only square matrices allowed
  {
    auto matrix = std::make_unique<Matrix<Scalar>>(dim, dim);
    auto vector = std::make_unique<Vector<Scalar>>(dim);

    MatrixVariant mv(std::move(*matrix));
    VectorVariant vv(std::move(*vector));

    return std::make_unique<MatrixEngine>(mv, vv);
  }

  //----------------------------------------------------------------------
  // Transition‐matrix postprocessing (used by assembleTransitionMatrix_)
  //----------------------------------------------------------------------
  void addIdentityInPlace()
  {
    std::visit(overloaded
    {
      [&](DoubleMatrixWrap &W){ W.addIdentity(); },
      [&](MPRealMatrixWrap &W){ W.addIdentity(); }
    }, matWrap_);
  }

  void pruneInPlace()
  {
    std::visit(overloaded
    {
      [&](DoubleMatrixWrap &W){ W.prune(); },
      [&](MPRealMatrixWrap &W){ W.prune(); }
    }, matWrap_);
  }

  void compressInPlace()
  {
    std::visit(overloaded
    {
      [&](DoubleMatrixWrap &W){ W.compress(); },
      [&](MPRealMatrixWrap &W){ W.compress(); }
    }, matWrap_);
  }

  //----------------------------------------------------------------------
  // Arithmetic & solve (wrapper‐level APIs)
  //----------------------------------------------------------------------
  MatrixEngine& operator+=(MatrixEngine const &rhs)
  {
    visitSameType(matWrap_, rhs.matWrap_, [&](auto &A, auto const &B)
    {
      A += B;
    });

    return *this;
  }

  friend MatrixEngine operator+(MatrixEngine lhs, MatrixEngine const &rhs)
  {
    lhs += rhs;
    return lhs;
  }

  MatrixEngine& operator*=(double s)
  {
    std::visit(overloaded
    {
      [&](DoubleMatrixWrap &M){ M *= s; },
      [&](MPRealMatrixWrap &M){ M *= mpfr::mpreal(s); }
    }, matWrap_);
    std::visit(overloaded
    {
      [&](DoubleVectorWrap &v){ v *= s; },
      [&](MPRealVectorWrap &v){ v *= mpfr::mpreal(s); }
    }, vecWrap_);
    return *this;
  }

  friend MatrixEngine operator*(MatrixEngine m, double s)
  {
    m *= s;
    return m;
  }

  friend MatrixEngine operator*(double s, MatrixEngine m)
  {
    m *= s;
    return m;
  }

  VectorVariant operator*(const VectorVariant &v) const
  {
    VectorVariant out;
    visitSameType(matWrap_, v,
                  [&](auto const &M, auto const &vec)
                  {
                    out = M * vec;
                  });
    return out;
  }

  // Solve A x = v at wrapper level
  VectorVariant solveSystem() const
  {
    VectorVariant sol;
    std::visit(overloaded
    {
      [&](DoubleMatrixWrap const &M){
        auto up = M.solve(std::get<DoubleVectorWrap>(vecWrap_));
        sol   = *up;
      },
      [&](MPRealMatrixWrap const &M){
        auto up = M.solve(std::get<MPRealVectorWrap>(vecWrap_));
        sol   = *up;
      }
    }, matWrap_);
    return sol;
  }

  /// scale only the transition matrix (not the vector)
  void scaleMatrix(double s)
  {
    std::visit(overloaded
    {
      [&](DoubleMatrixWrap &M){ M *= s; },
      [&](MPRealMatrixWrap &M){ M *= mpfr::mpreal(s); }
    }, matWrap_);
  }

  /// in‐place add another raw MatrixVariant
  void addToMatrix(const MatrixVariant &other)
  {
    visitSameType(matWrap_, other, [&](auto &M, auto const &N)
    {
      M += N;
    });
  }

  // prints type held by instance of MatrixVariant
  std::string getMatrixType() const
  {
    switch(matWrap_.index())
    {
     case 0: return "double";
     case 1: return "mpreal";
     default: throw bpp::Exception("MatrixEngine::unexpected MatrixVariant index");
    }
  }

  // prints type held by instance of VectorVariant
  std::string getVectorType() const
  {
    switch(vecWrap_.index())
    {
     case 0: return "double";
     case 1: return "mpreal";
     default: throw bpp::Exception("MatrixEngine::unexpected VectorVariant index");
    }
  }

private:
  //----------------------------------------------------------------------
  // Overwrite wrappers from raw Eigen variants (used by Epoch after solve)
  //----------------------------------------------------------------------
  void setMatrixFromEigen(MatrixEigenVariant newM)
  {
    MatrixVariant converted = std::visit(overloaded
    {
      [](DoubleMatrixEigen const &Me) -> MatrixVariant
      {
        return MatrixVariant{ DoubleMatrixWrap(Me) };
      },

      [](MPRealMatrixEigen const &Me) -> MatrixVariant
      {
        return MatrixVariant{ MPRealMatrixWrap(Me) };
      }
    }, std::move(newM));

    matWrap_ = std::move(converted);
    std::visit(overloaded{[&](auto const &M)
    {
      rows_ = M.rows();
      cols_ = M.cols();
    }
    }, matWrap_);
  }

  void setVectorFromEigen(VectorEigenVariant newV)
  {
    VectorVariant converted = std::visit(overloaded{
      [](DoubleVectorEigen const &ve) -> VectorVariant
      {
        return VectorVariant{ DoubleVectorWrap(ve) };
      },

      [](MPRealVectorEigen const &ve) -> VectorVariant
      {
        return VectorVariant{ MPRealVectorWrap(ve) };
      }
    }, std::move(newV));

    vecWrap_ = std::move(converted);
    int sz = std::visit(overloaded
    {
      [&](DoubleVectorWrap const &V){ return V.size(); },
      [&](MPRealVectorWrap const &V){ return V.size(); }
    }, vecWrap_);

    if (sz != cols_)
      throw bpp::Exception("MatrixEngine::setVectorFromEigen: size mismatch");
  }
};

#endif // _MATRIXENGINE_HPP_
