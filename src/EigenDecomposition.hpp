/*
 * Authors: Gustavo V. Barroso
 * Created: 06/09/2022
 * Last modified: 12/09/2022
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
  Eigen::EigenSolver<Eigen::MatrixXd> es_;
  Eigen::MatrixXd matInverse_;
  Eigen::MatrixXd lambda_;
  Eigen::MatrixXd mat_;

public:
  EigenDecomposition():
  es_(),
  matInverse_(),
  lambda_(),
  mat_()
  { }

  void exponentiate(Eigen::MatrixXd& mat, size_t exponent)
  {
    es_.compute(mat);

    matInverse_ = es_.eigenvectors().real().inverse();
    lambda_ = es_.eigenvalues().real().array().pow(exponent).matrix().asDiagonal(); // NOTE
    mat_ = es_.eigenvectors().real();

    mat = matInverse_ * lambda_ * mat_;
  }

  const Eigen::MatrixXd& matrix()
  {
    return mat_;
  }

  const Eigen::MatrixXd& matrixInverse()
  {
    return matInverse_;
  }

  const Eigen::MatrixXd& lambdaReal()
  {
    return lambda_;
  }

  void setLambda(const Eigen::MatrixXd& newLambda)
  {
    lambda_ = newLambda;
  }

};

#endif
