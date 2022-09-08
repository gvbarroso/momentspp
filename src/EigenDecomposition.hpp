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
  Eigen::MatrixXd matInverse_;
  Eigen::MatrixXd lambda_;
  Eigen::MatrixXd mat_;

public:
  EigenDecomposition():
  mat_(),
  lambda_(),
  matInverse_()
  { }

  EigenDecomposition(const Eigen::MatrixXd& mat, size_t exponent):
  mat_(),
  lambda_(),
  matInverse_()
  {
    Eigen::EigenSolver<Eigen::MatrixXd> es(mat);

    matInverse_ = es.eigenvectors().inverse();
    lambda_ = es.eigenvalues().pow(exponent); // NOTE
    mat_ = es.eigenvectors();
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

  void setLambda(const Eigen::MatrixXd& newLambda)
  {
    lambda_ = newLambda;
  }

};

#endif
