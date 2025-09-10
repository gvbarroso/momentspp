/*
 * Authors: Gustavo V. Barroso
 * Created: 08/09/2025
 * Last modified: 10/09/2025
 *
 */

#ifndef _MATDOUBLEH_
#define _MATDOUBLEH_

#include "MatrixInterface.hpp"

class MatrixDouble: public MatrixInterface
{
private:
  Eigen::SparseMatrix<double> mat_;

public:
  MatrixDouble(size_t rows, size_t cols):
  mat_(rows, cols)
  { }

  MatrixDouble(const Eigen::SparseMatrix<double>& mat):
  mat_(mat)
  { }

  MatrixDouble(Eigen::SparseMatrix<double>&& mat):
  mat_(std::move(mat))
  { }

  const Eigen::SparseMatrix<double>& const getMatrix()
  {
    return mat_;
  }

  void setMatrix(const Eigen::SparseMatrix<double>& mat)
  {
    mat_ = mat;
  }

  void insert(size_t row, size_t col, double value) override
  {
    mat_.insert(row, col) = value;
  }

  void scale(double scalar) override
  {
    for(size_t i = 0; i < mat.outerSize(); ++i)
    {
      for(Eigen::SparseMatrix<double>::InnerIterator it(mat, i); it; ++it)
        it.valueRef() *= scalar;
    }
  }

  void zeroNegatives() override
  {
    for(size_t i = 0; i < mat_.outerSize(); ++i)
    {
      for(Eigen::SparseMatrix<double>::InnerIterator it(mat_, i); it; ++it)
      {
        if(it.value() < 0)
          it.valueRef() = 0.0;
      }
    }
    mat_.prune(0.0);
  }

  void print(const std::string& fileName) const override
  {
    std::ofstream matFile;
    matFile.open(fileName);

    for(size_t i = 0; i < mat_.rows(); ++i)
    {
      for(size_t j = 0; j < mat_.cols(); ++j)
      {
        matFile << mat_.coeffRef(i, j);

        if(j < mat_.cols() - 1)
          matFile << ",";
      }

      matFile  << "\n";
    }

    matFile.close();
  }

  std::unique_ptr<MatrixInterface> clone()
  {
    return std::make_unique<MatrixDouble>(mat_);
  }

  void addInPlace(const MatrixInterface& other) override
  {
    const auto& otherMat = static_cast<const MatrixDouble&>(other).getMatrix();

    if(mat_.rows() != otherMat.rows() || mat_.cols() != otherMat.cols())
      throw bpp::Exception("matrix dimensions do not match!");

    mat_ += otherMat;
    mat_.makeCompressed();
  }

  std::unique_ptr<VectorInterface> multiply(const VectorInterface& vec) const override
  {
    const Eigen::VectorXd& v = static_cast<const VectorDouble&>(vec).data();
    Eigen::VectorXd result = mat_ * v;

    auto out = std::make_unique<VectorDouble>(result.size());
    for(size_t i = 0; i < result.size(); ++i)
      out->set(i, result(i));

    return out;
  }

  std::unique_ptr<MatrixInterface> multiply(const MatrixInterface& other) const override
  {
    const auto& otherMat = static_cast<const MatrixDouble&>(other).getMatrix();

    Eigen::SparseMatrix<double> result = mat_ * otherMat;
    return std::make_unique<MatrixDouble>(std::move(result));
  }

  std::unique_ptr<MatrixInterface> add(const MatrixInterface& other) const override
  {
    const auto& otherMat = static_cast<const MatrixDouble&>(other.getMatrix();

    Eigen::SparseMatrix<double> result = mat_ + otherMat;
    return std::make_unique<MatrixDouble>(std::move(result));
  }

  void setFromTriplets(const std::vector<Eigen::Triplet<double>>& triplets) override
  {
    mat_.setFromTriplets(triplets.begin(), triplets.end());
  }

  void makeCompressed() override
  {
    mat_.makeCompressed();
  }

  std::unique_ptr<MatrixInterface> identity() const override
  {
    size_t n = mat_.rows();
    auto I = std::make_unique<MatrixDouble>(n, n);

    std::vector<Eigen::Triplet<double>> triplets;
    for(size_t i = 0; i < n; ++i)
      triplets.emplace_back(i, i, 1.0);

    I->mat_.setFromTriplets(triplets.begin(), triplets.end());
    I->mat_.makeCompressed();

    return I;
  }

  std::unique_ptr<VectorInterface> solve(const VectorInterface& rhs) const override
  {
    const auto& b = static_cast<const VectorDouble&>(rhs).data();

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(mat);

    if(solver.info() != Eigen::Success)
      throw bpp::Exception("Matrix decomposition failed!\n");

    Eigen::VectorXd x = solver.solve(b);

    if(solver.info() != Eigen::Success)
      throw bpp::Exception("Solving linear system failed!\n");

    auto result = std::make_unique<VectorDouble>(x.size());
    for(size_t i = 0; i < x.size(); ++i)
      result->set(i, x(i));

    return result;
  }

};

#endif
