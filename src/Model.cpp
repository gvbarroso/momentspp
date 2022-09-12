/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 12/09/2022
 *
 */


#include "Model.hpp"

void Model::fireParameterChanged(const bpp::ParameterList& params)
{
  matchParametersValues(params);
  updateEpochs_(params); // updates transitionMatrix_ within each epoch

  computeExpectedSumStats_();
  computeCompositeLogLikelihood_(sslib_.getYvec(), sslib_.getCovarMatrix()); // for each rec. bin
}

void Model::updateEpochs_(const bpp::ParameterList& params)
{
  for(auto it = std::begin(epochs_); it != std::end(epochs_); ++it)
    (*it)->fireParameterChanged(params);
}

void Model::computeExpectedSumStats_()
{
  expected_ = steadYstate_; // resets "to the deep past"

  // TODO use recipe of pop splits and admixs to change vector of sum stats between each epoch
  for(auto it = std::begin(epochs_); it != std::end(epochs_); ++it) // epochs must be sorted from past to present
    (*it)->computeExpectedSumStats(expected_); // trickling sum stats vector down the epochs
}

void Model::computeSteadyState()
{
  // we find the eigenvector associated with lambda == 1 in the (arbitrary) transition matrix from epoch 0
  Eigen::EigenSolver<Eigen::MatrixXd> es(epochs_[0]->getTransitionMatrix());
  steadYstate_ = es.eigenvectors().real().col(0);

  std::cout << "leading eigenvalue = " << es.eigenvalues()[0] << "; associated eigenvalues: " << steadYstate_ << std::endl;
}

void Model::computeCompositeLogLikelihood_(const Eigen::VectorXd& obsMeans, const Eigen::MatrixXd& obsCovarMat)
{
  double cll = 0.;

  /*for(auto it = std::begin(recBins_); it != std::end(recBins_); ++it)
  {
    cll += det(2*covarMat)^(-1/2) * exp(-1/2 * (expected_ - means).transpose() * covarMat^(-1)*(expected_ - means);
  }*/

  compLogLikelihood_ = cll;
}
