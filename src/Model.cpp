/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 02/09/2022
 *
 */


#include "Model.hpp"

void Model::fireParameterChanged(const bpp::ParameterList& params)
{
  matchParametersValues(params);
  updateEpochs_(params);

  auto expectedStats = computeExpectedSumStats_();

  computeCompositeLogLikelihood_(sslib_->getYvec(), expectedStats);
}

void Model::updateEpochs_(const bpp::ParameterList& params)
{
  for(auto it = std::begin(epochs_); it != std::end(epochs_); ++it)
    (*it)->fireParametersChanged(params);
}

void Model::computeCompositeLogLikelihood_(const Eigen::Matrix<double, Dynamic, 1>& observed, const Eigen::Matrix<double, Dynamic, 1>& expected)
{
  // compLogLikelihood_ = f(observed, expected) for different windows or something
}

Eigen::Matrix<double, Dynamic, 1> Model::computeExpectedSumStats_()
{
  // NOTE mind lexicographic order in SumStatsLibrary
  Eigen::Matrix<double, Dynamic, 1> y = fetchSteadyState();

  for(auto it = std::begin(epochs_); it != std::end(epochs_); ++it) // epochs must be sorted from past to present
    it->computeExpectedSumStats(y); // trickling sum stats vector down the epochs

  return y;
}

Eigen::Matrix<double, Dynamic, 1> Model::fetchSteadyState()
{
    // find eigenvector associated with lambda == 1
}


