/*
 * Authors: Gustavo V. Barroso
 * Created: 31/08/2022
 * Last modified: 07/09/2022
 *
 */


#include "Epoch.hpp"

void Epoch::fireParameterChanged(const bpp::ParameterList& params)
{
  if(matchParametersValues(params))
    updateOperators_(params);
}

void Epoch::updateOperators_(const bpp::ParameterList& params)
{
  for(auto it = std::begin(operators_); it != std::end(operators_); ++it)
    (*it)->fireParametersChanged(params, duration());
}

Eigen::MatrixXd Epoch::fetchCombinedOperators()
{
  Eigen::MatrixXd mat(sslib_.getNumStats(), sslib_.getNumStats());
  mat.setIdentity();

  // NOTE we must be careful with the order of operations
  for(auto it = std::begin(operators_); it != std::end(operators_); ++it)
    mat *= (*it)->fetchCombinedMatrix(duration());

  return mat;
}

void Epoch::computeExpectedSumStats(Eigen::VectorXd& y)
{
  combineOperators() * y;
}
