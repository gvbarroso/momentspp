/*
 * Authors: Gustavo V. Barroso
 * Created: 31/08/2022
 * Last modified: 23/09/2022
 *
 */


#include "Epoch.hpp"

// this method is where the heavier Eigen linear algebra takes place
void Epoch::fireParameterChanged(const bpp::ParameterList& params)
{
  if(matchParametersValues(params))
    updateOperators_(params);

  Eigen::SparseMatrix<double> mat = operators_[0]->fetchCombinedMatrix(); // init

  // NOTE we must be careful with the order of operations
  for(size_t i = 1; i < operators_.size(); ++i)
    mat = mat * operators_[i]->fetchCombinedMatrix();

  transitionMatrix_ = mat; // converts to dense format
  eigenDec_.exponentiate(transitionMatrix_, duration()); // matrix passed as non-const ref
}

void Epoch::computeSteadyState_()
{
  updateOperators_(getParameters());
  Eigen::SparseMatrix<double> mat = operators_[0]->fetchCombinedMatrix(); // init

  // NOTE we must be careful with the order of operations
  for(size_t i = 1; i < operators_.size(); ++i)
    mat = mat * operators_[i]->fetchCombinedMatrix();

  transitionMatrix_ = mat; // converts to dense format

  // we find the eigenvector associated with (leading) eigenvalue == 1 in transitionMatrix_
  Eigen::EigenSolver<Eigen::MatrixXd> es(transitionMatrix_);
  steadYstate_ = es.eigenvectors().real().col(0);

  std::cout << "leading eigenvalue = " << es.eigenvalues()[0] << "; associated eigenvector: " << steadYstate_ << std::endl;
}

void Epoch::transferStatistics(Eigen::VectorXd& y)
{
  Eigen::VectorXd tmp(ssl_.getMoments().size()); // y and tmp have potentially different sizes

  // for each Moment in *this epoch, we assign its value from its parental Moment from the previous epoch, from which y comes
  for(size_t i = 0; i < tmp.size(); ++i)
    tmp(i) = y[ssl_.getMoments()[i].getParent()->getPosition()];

  y = tmp; // swap
}

void Epoch::updateMoments(const Eigen::VectorXd& y)
{
  if(y.size() != ssl_.getMoments().size())
    throw bpp::Exception("Epoch::attempted to update moments from vector of different size!");

  for(size_t i = 0; i < y.size(); ++i)
    ssl_.getMoments()[i].setValue(y(i));
}
