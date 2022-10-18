/*
 * Authors: Gustavo V. Barroso
 * Created: 31/08/2022
 * Last modified: 23/09/2022
 *
 */

#include <ios>

#include "Epoch.hpp"

// this method is where the heavier Eigen linear algebra takes place
void Epoch::fireParameterChanged(const bpp::ParameterList& params)
{
  if(matchParametersValues(params))
    updateOperators_(params);

  Eigen::SparseMatrix<double> mat = operators_[0]->fetchCombinedMatrix(); // init mat

  // we must be careful with the order of operations
  for(size_t i = 1; i < operators_.size(); ++i)
    mat = mat * operators_[i]->fetchCombinedMatrix();

  transitionMatrix_ = mat; // converts to dense format
  eigenDec_.exponentiate(transitionMatrix_, duration()); // matrix passed as non-const ref
}

void Epoch::computeSteadyState_()
{
  //updateOperators_(getParameters());
  Eigen::SparseMatrix<double> mat = operators_[0]->fetchCombinedMatrix(); // init mat
  //operators_[0]->getParameters().printParameters(std::cout);
  std::cout << std::scientific << mat << std::endl;

  // we must be careful with the order of operations
  for(size_t i = 1; i < operators_.size(); ++i)
  {
    auto m = operators_[i]->fetchCombinedMatrix();
    mat = m * mat;
    //operators_[i]->getParameters().printParameters(std::cout);
    std::cout << std::scientific << m << std::endl;
    std::cout << std::scientific << mat << std::endl;
  }

  transitionMatrix_ = mat; // converts to dense format

  // we find the eigenvector associated with (leading) eigenvalue == 1 in transitionMatrix_
  Eigen::EigenSolver<Eigen::MatrixXd> es(transitionMatrix_);
  for(int i = 0; i < es.eigenvalues().size(); ++i)
  {
    std::cout << "\n" << std::setprecision(32) << es.eigenvalues().real()(i) << "\n";
    std::cout << es.eigenvectors().col(i) << "\n\n";

    if(es.eigenvalues().real()(i) == 1.)
      steadYstate_ = es.eigenvectors().col(i).real();
  }

  //std::cout << std::setprecision(7) << transitionMatrix_ << "\n";
  std::cout << std::setprecision(7) << steadYstate_ << "\n";
  //std::cout << es.eigenvectors().real() << "\n";
  //std::cout << es.eigenvalues().real() << "\n";
}

void Epoch::transferStatistics(Eigen::VectorXd& y)
{
  Eigen::VectorXd tmp(ssl_.getMoments().size()); // y and tmp have potentially different sizes
  tmp.setZero();

  // for each Moment in *this epoch, we assign its value from its parental Moment from the previous epoch (from which y comes)
  for(size_t i = 0; i < tmp.size(); ++i)
  {
    int idx = ssl_.getMoments()[i].getParent()->getPosition();
    tmp(i) = y(idx);
  }

  y = tmp; // swap
}

void Epoch::updateMoments(const Eigen::VectorXd& y)
{
  if(y.size() != ssl_.getMoments().size())
    throw bpp::Exception("Epoch::attempted to update moments from vector of different size!");

  for(size_t i = 0; i < y.size(); ++i)
    ssl_.getMoments()[i].setValue(y(i));
}
