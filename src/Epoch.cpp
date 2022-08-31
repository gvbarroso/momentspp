/*
 * Authors: Gustavo V. Barroso
 * Created: 31/08/2022
 * Last modified: 31/08/2022
 *
 */


#include "Epoch.hpp"

void Epoch::fireParameterChanged(const bpp::ParameterList& params, Eigen::Matrix<double, Dynamic, 1>& y)
{
  if(matchParametersValues(params))
    updateOperators_(params);
}

void Epoch::updateOperators_(const bpp::ParameterList& params)
{
  for(auto it = std::begin(operators_); it != std::end(operators_); ++it)
    (*it)->fireParametersChanged(params);
}

Eigen::Matrix<double, Dynamic, Dynamic> Epoch::integrateOperators_()
{
  Eigen::Matrix<double, Dynamic, Dynamic> matrix(sslib_.getNumStats(), sslib_.getNumStats());

  if(continuousTime_) // we combine operators by matrix addition
    for(auto it = std::begin(operators_); it != std::end(operators_); ++it)
      matrix += (*it)->fetchCombinedMatrix(duration());

  else // we combine operators by matrix multiplication WARNING we must be careful with the order of operations
    for(auto it = std::begin(operators_); it != std::end(operators_); ++it)
      matrix *= (*it)->fetchCombinedMatrix(duration());

  return matrix;
}

void Epoch::computeExpectedSumStats(Eigen::Matrix<double, Dynamic, 1>& y)
{
  integrateOperators_() * y;
}
