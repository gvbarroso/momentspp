/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 05/09/2022
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
  auto y = steadYstate_;

  for(auto it = std::begin(epochs_); it != std::end(epochs_); ++it) // epochs must be sorted from past to present
    it->computeExpectedSumStats(y); // trickling sum stats vector down the epochs

  return y;
}

void Model::computeSteadyState()
{
    // we find the eigenvector associated with lambda == 1 in an arbitrary combined matrix
    auto matrix = epochs_[0]->fetchCombinedOperators();
    Eigen::EigenSolver es(matrix);

    steadYstate_ = es.eigenvectors().col(0);

    std::cout << "leading eigenvalue = " << es.eigenvalues()[0] << "; associated eigenvalues: " << steadYstate_ << std::endl;

}


