/*
 * Authors: Gustavo V. Barroso
 * Created: 31/08/2022
 * Last modified: 11/09/2022
 *
 */


#include "Epoch.hpp"

void Epoch::fireParameterChanged(const bpp::ParameterList& params)
{
  if(matchParametersValues(params))
    updateOperators_(params);

  Eigen::SparseMatrix<double> mat(operators_[0]->getMatrix()[0].rows(), operators_[0]->getMatrix()[0].cols());
  for(size_t i = 0; i < mat.cols(); ++i) // mat is square
    mat(i, i) += 1.;

  // NOTE we must be careful with the order of operations
  for(auto it = std::begin(operators_); it != std::end(operators_); ++it)
    mat *= (*it)->fetchCombinedMatrix();

  transitionMatrix_.setIdentity();
  transitionMatrix_ *= mat;

  eigenDec_.exponentiate(transitionMatrix_, duration());
}
