/*
 * Authors: Gustavo V. Barroso
 * Created: 06/09/2022
 * Last modified: 08/09/2022
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
  Eigen::MatrixXd matInverse_;
  Eigen::MatrixXd lambda_;
  Eigen::MatrixXd mat_;

public:
  EigenDecomposition():
  matInverse_(),
  lambda_(),
  mat_()
  { }

  EigenDecomposition(const Eigen::MatrixXd& mat, size_t exponent):
  matInverse_(),
  lambda_(),
  mat_()
  {
    Eigen::EigenSolver<Eigen::MatrixXd> es(mat);

    matInverse_ = es.eigenvectors().real().inverse();
    lambda_ = es.eigenvalues().real().array().pow(exponent).matrix().asDiagonal(); // NOTE
    mat_ = es.eigenvectors().real();
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
