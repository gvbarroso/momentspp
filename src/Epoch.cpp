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

  Eigen::SparseMatrix<double> mat = operators_[0]->fetchCombinedMatrix(); // init

  // NOTE we must be careful with the order of operations
  for(size_t i = 1; i < operators_.size(); ++i)
    mat *= operators_[i]->fetchCombinedMatrix();

  transitionMatrix_ = mat; // converts to dense format
  eigenDec_.exponentiate(transitionMatrix_, duration());
}
