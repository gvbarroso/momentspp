/*
 * Authors: Gustavo V. Barroso
 * Created: 08/09/2025
 * Last modified: 10/09/2025
 *
 */

#ifndef _MATMPREAL_H_
#define _MATMPREAL_H_

#include "MatrixInterface.hpp"

class MatrixMPReal: public MatrixInterface
{
private:
  Eigen::SparseMatrix<mpfr::mpreal> mat_;

public:
  MatrixMPReal(size_t rows, size_t cols):
  mat_(rows, cols)
  { }

  MatrixMPReal(const Eigen::SparseMatrix<mpfr::mpreal>& mat):
  mat_(mat)
  { }

  MatrixMPReal(Eigen::SparseMatrix<mpfr::mpreal>&& mat):
  mat_(std::move(mat))
  { }

  const Eigen::SparseMatrix<mpfr::mpreal>& const getMatrix()
  {
    return mat_;
  }

  void setMatrix(const Eigen::SparseMatrix<mpfr::mpreal>& mat)
  {
    mat_ = mat;
  }

  void insert(size_t row, size_t col, double value) override
  {
    mat_.insert(row, col) = value;
  }

  void scale(double scalar) override
  {
    mpfr::mpreal s = mpfr::mpreal(scalar);
    for(size_t i = 0; i < mat.outerSize(); ++i)
    {
      for(Eigen::SparseMatrix<mpfr::mpreal>::InnerIterator it(mat, i); it; ++it)
        it.valueRef() *= s;
    }
  }

  void scale(const mpfr::mpreal& scalar)
  {
    for(size_t i = 0; i < mat.outerSize(); ++i)
    {
      for(Eigen::SparseMatrix<mpfr::mpreal>::InnerIterator it(mat, i); it; ++it)
        it.valueRef() *= scalar;
    }
  }

  void print(const std::string& fileName) const override
  {
    std::ofstream matFile;
    matFile.open(fileName);

    for(size_t i = 0; i < mat_.rows(); ++i)
    {
      for(size_t j = 0; j < mat_.cols(); ++j)
      {
        mpfr::mpreal value = mat.coeff(i, j);
        matFile << value.toString();

        if(j < mat_.cols() - 1)
          matFile << ",";
      }

      matFile  << "\n";
    }

    matFile.close();
  }

  std::unique_ptr<MatrixInterface> clone()
  {
    return std::make_unique<MatrixMPReal>(mat_);
  }

  void addInPlace(const MatrixInterface& other) override
  {
    const auto& otherMat = static_cast<const MatrixMPReal&>(other).getMatrix();

    if(mat_.rows() != otherMat.rows() || mat_.cols() != otherMat.cols())
      throw bpp::Exception("matrix dimensions do not match!");

    mat_ += otherMat;
    mat_.makeCompressed();
  }

  std::unique_ptr<VectorInterface> multiply(const VectorInterface& vec) const override
  {
    const auto& v = static_cast<const VectorMPReal&>(vec).data();
    auto result = mat * v;

    auto out = std::make_unique<VectorMPReal>(result.size());
    for(size_t i = 0; i < result.size(); ++i)
      out->set(i, result(i));

    return out;
  }

  std::unique_ptr<MatrixInterface> multiply(const MatrixInterface& other) const override
  {
    const MatrixMPReal* otherMPReal = dynamic_cast<const MatrixMPReal*>(&other);
    if(!otherMPReal) throw bpp::Exception("type mismatch in multiplication");

    Eigen::SparseMatrix<mpfr::mpreal> result = mat_ * otherMPReal->getMatrix();
    auto product = std::make_unique<MatrixMPReal>(result.rows(), result.cols());
    product->setMatrix(result);

    return product;
  }

  std::unique_ptr<MatrixInterface> add(const MatrixInterface& other) const override
  {
    const MatrixMPReal* otherMPReal = dynamic_cast<const MatrixMPReal*>(&other);
    if(!otherMPReal) throw bpp::Exception("type mismatch in addition!");

    Eigen::SparseMatrix<mpfr::mpreal> result = mat_ + otherMPReal->getMatrix();
    auto sum = std::make_unique<MatrixMPReal>(result.rows(), result.cols());
    sum->setMatrix(result);

    return sum;
  }

  void setFromTriplets(const std::vector<Eigen::Triplet<double>>& triplets) override
  {
    std::vector<Eigen::Triplet<mpfr::mpreal>> mpTriplets;
    mpTriplets.reserve(triplets.size());

    for(const auto& t : triplets)
      mpTriplets.emplace_back(t.row(), t.col(), mpfr::mpreal(t.value()));

    mat_.setFromTriplets(mpTriplets.begin(), mpTriplets.end());
  }

  void makeCompressed() override
  {
    mat_.makeCompressed();
  }

  std::unique_ptr<VectorInterface> solve(const VectorInterface& rhs) const override
  {
    const auto& b = static_cast<const VectorMPReal&>(rhs).data();

    // decomposing a MPReal matrix requires dense format
    Eigen::Matrix<mpfr::mpreal, Eigen::Dynamic, Eigen::Dynamic> dense = mat_.toDense();
    Eigen::FullPivLU<Eigen::Matrix<mpfr::mpreal, Eigen::Dynamic, Eigen::Dynamic>> solver(dense);

    if(!solver.isInvertible())
      throw bpp::Exception("Matrix is not invertible!\n");


    Eigen::Matrix<mpfr::mpreal, Eigen::Dynamic, 1> x = solver.solve(b);

    auto result = std::make_unique<VectorMPReal>(x.size());
    for(size_t i = 0; i < x.size(); ++i)
        result->set(i, x(i));

    return result;
  }

  Eigen::Matrix<mpfr::mpreal, Eigen::Dynamic, Eigen::Dynamic> toDense() const
  {
    Eigen::Matrix<mpfr::mpreal, Eigen::Dynamic, Eigen::Dynamic> dense(mat_.rows(), mat_.cols());
    dense.setZero();

    for(size_t i = 0; i < mat_.outerSize(); ++i)
      for(Eigen::SparseMatrix<mpfr::mpreal>::InnerIterator it(mat_, i); it; ++it)
        dense(it.row(), it.col()) = it.value();

    return dense;
  }

};

#endif
