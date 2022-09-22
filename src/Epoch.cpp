/*
 * Authors: Gustavo V. Barroso
 * Created: 31/08/2022
 * Last modified: 22/09/2022
 *
 */


#include "Epoch.hpp"

// NOTE this method is where the heavier linear algebra is
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
