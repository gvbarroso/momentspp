/*
 * Authors: Gustavo V. Barroso
 * Created: 06/09/2022
 * Last modified: 07/11/2022
 *
 */


#ifndef _EIGENDECOMP_H_
#define _EIGENDECOMP_H_


#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "Log.hpp"

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

  void exponentiate(Eigen::MatrixXd& mat, size_t exponent)
  {
    Log logger;
    logger.openFile("input_mat.txt");
    logger.getLogFile() << mat << "\n\n";
    logger.closeFile();

    es_.compute(mat);

    logger.openFile("dec_mat.txt");
    logger.getLogFile() << es_.eigenvectors().inverse() * es_.eigenvalues().matrix().asDiagonal() * es_.eigenvectors() << "\n\n";
    logger.closeFile();

    // WARNING
    mat_ = es_.eigenvectors().real();
    matInverse_ = es_.eigenvectors().real().inverse();
    lambda_ = es_.eigenvalues().real().array().pow(exponent).matrix().asDiagonal();

    std::cout << "\n\neigenvalues:\n" << es_.eigenvalues();

    mat = mat_ * lambda_ * matInverse_;

    logger.openFile("dec_mat_real.txt");
    logger.getLogFile() << matInverse_ << "\n\n";
    logger.getLogFile() << lambda_ << "\n\n";
    logger.getLogFile() << mat_ << "\n\n";
    logger.getLogFile() << matInverse_ * mat_  << "\n\n";
    logger.getLogFile() << lambda_ * mat_ * matInverse_ << "\n\n";
    logger.closeFile();
  }

};

#endif
