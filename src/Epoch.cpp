/*
 * Authors: Gustavo V. Barroso
 * Created: 31/08/2022
 * Last modified: 08/09/2022
 *
 */


#include "Epoch.hpp"

void Epoch::fireParameterChanged(const bpp::ParameterList& params)
{
  if(matchParametersValues(params))
    updateOperators_(params);

  transitionMatrix_.setIdentity();

  // NOTE we must be careful with the order of operations
  for(auto it = std::begin(operators_); it != std::end(operators_); ++it)
    transitionMatrix_ *= (*it)->fetchCombinedMatrix();
}

void Epoch::updateOperators_(const bpp::ParameterList& params)
{
  for(auto it = std::begin(operators_); it != std::end(operators_); ++it)
    (*it)->fireParameterChanged(params);
}

void Epoch::computeExpectedSumStats(Eigen::VectorXd& y)
{
  transitionMatrix_ * y;
}
