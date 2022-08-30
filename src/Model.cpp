/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 25/08/2022
 *
 */


#include "Model.hpp"

void Model::fireParameterChanged(const bpp::ParameterList& params)
{
  matchParametersValues(params);
  updateOperators_(params);

  auto combinedOperator = integrateOperators_();
  auto expectedStats = computeExpectedSumStats_(combinedOperator);

  computeCompositeLogLikelihood_(sslib_->getYvec(), expectedStats);
  computeAic_();
}

void Model::updateOperators_(const bpp::ParameterList& params)
{
  for(auto it = std::begin(operators_); it != std::end(operators_); ++it)
    (*it)->fireParametersChanged(params);
}

Eigen::Matrix<double, Dynamic, Dynamic> Model::integrateOperators_()
{
  Eigen::Matrix<double, Dynamic, Dynamic> matrix(sslib_.getNumStats(), sslib_.getNumStats());

  if(continuousTime_) // we combine operators by matrix addition
    for(auto it = std::begin(operators_); it != std::end(operators_); ++it)
      matrix += (*it)->getCombinedMatrix();

  else // we combine operators by matrix multiplication WARNING we must be careful with the order of operations
    for(auto it = std::begin(operators_); it != std::end(operators_); ++it)
      matrix *= (*it)->getCombinedMatrix();

  return matrix;
}

void Model::computeCompositeLogLikelihood_(const Eigen::Matrix<double, Dynamic, 1>& observed, const Eigen::Matrix<double, Dynamic, 1>& expected)
{
  // logLikelihood_ = f(observed, expected) for different windows or something
}

Eigen::Matrix<double, Dynamic, 1> Model::computeExpectedSumStats_(const Eigen::Matrix<double, Dynamic, Dynamic>& matrix)
{
  Eigen::EigenSolver es(matrix); // NOTE put one es inside each epoch/operator?

  // we want to do something like this:
  Eigen::Matrix<double, Dynamic, 1> expected = (es.eigenvectors() * es.eigenvalues() ^ t * es.eigenvectors().inverse()) * data_;

  return expected;
}

