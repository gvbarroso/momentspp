/*
 * Authors: Gustavo V. Barroso
 * Created: 06/09/2022
 * Last modified: 06/09/2022
 *
 */


#ifndef _EIGENDECOMP_H_
#define _EIGENDECOMP_H_


#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>


class EigenDecomposition:
{

private:
  Eigen::EigenSolver<MatrixType>::eigenvectors mat_;
  Eigen::EigenSolver<MatrixType_>::eigenvalues lambda_;
  Eigen::EigenSolver<MatrixType>::eigenvectors matInverse_;

public:
  EigenDecomposition(const Eigen::EigenSolver& es):
  mat_(es.eigenvectors()),
  lambda_(es.eigenvalues()),
  matInverse_(es.eigenvectors().inverse())
  { }

  EigenDecomposition(const Eigen::SparseMatrix<double> mat& mat):
  mat_(),
  lambda_(),
  matInverse_()
  {
    Eigen::EigenSolver es(mat);
    mat_ = es.eigenvectors();
    lambda_ = es.eigenvalues();
    matInverse_ = es.eigenvectors().inverse();
  }

  const Eigen::Matrix& matrix()
  {
    return mat_;
  }

  const Eigen::Matrix& matrixInverse()
  {
    return matInverse_;
  }

  const Eigen::Matrix& lambda()
  {
    return lambda_.asDiagonal();
  }

};
