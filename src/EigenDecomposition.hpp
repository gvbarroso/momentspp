/*
 * Authors: Gustavo V. Barroso
 * Created: 06/09/2022
 * Last modified: 07/09/2022
 *
 */


#ifndef _EIGENDECOMP_H_
#define _EIGENDECOMP_H_


#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>


class EigenDecomposition
{

private:
  Eigen::MatrixXd mat_;
  Eigen::MatrixXd lambda_;
  Eigen::MatrixXd matInverse_;

public:
  EigenDecomposition(const Eigen::EigenSolver<Eigen::MatrixXd>& es):
  mat_(es.eigenvectors()),
  lambda_(es.eigenvalues()),
  matInverse_(es.eigenvectors().inverse())
  { }

  EigenDecomposition(const Eigen::MatrixXd& mat):
  mat_(),
  lambda_(),
  matInverse_()
  {
    Eigen::EigenSolver<Eigen::MatrixXd> es(mat);
    mat_ = es.eigenvectors();
    lambda_ = es.eigenvalues();
    matInverse_ = es.eigenvectors().inverse();
  }

  const Eigen::MatrixXd& matrix()
  {
    return mat_;
  }

  const Eigen::MatrixXd& matrixInverse()
  {
    return matInverse_;
  }

  const Eigen::MatrixXd& lambda()
  {
    return lambda_;
  }

  Eigen::MatrixXd lambdaMat(size_t exponent = 1)
  {
    return (lambda_.pow(exponent)).asDiagonal();
  }

  void setLambda(const Eigen::MatrixXd& newLambda)
  {
    lambda_ = newLambda;
  }

};

#endif
