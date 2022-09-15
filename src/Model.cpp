/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 14/09/2022
 *
 */


#include "Model.hpp"

void Model::fireParameterChanged(const bpp::ParameterList& params)
{
  matchParametersValues(params);
  updateEpochs_(params); // updates transitionMatrix_ within each epoch

  computeExpectedSumStats_();
  computeCompositeLogLikelihood_(data_.getYvec(), data_.getCovarMatrix()); // for each rec. bin
}

void Model::updateEpochs_(const bpp::ParameterList& params)
{
  for(auto it = std::begin(epochs_); it != std::end(epochs_); ++it)
    (*it)->fireParameterChanged(params);
}

void Model::computeExpectedSumStats_()
{
  expected_ = epochs_[0]->getSteadyState(); // resets stats to the "deep past"

  // we use map of pop (splits and admixtures) to change vector of sum stats between each epoch
  for(size_t i = 0; i < epochs_.size() - 1; ++i) // epochs must be sorted from past to present
  {
    epochs_[i]->computeExpectedSumStats(expected_); // trickling sum stats vector down the epochs

    // compare population indices between adjacent epochs

    auto from = it->getPops();
    auto to = (it)->getPops();
  }

  // final epoch
  epochs_.back()->computeExpectedSumStats(expected_);
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

void Model::popSplit_(const std::pair<size_t, std::pair<size_t, size_t>>& popTrio)
{

}

// admixes p2 and p3 (from second)
void Model::popAdmix_(const std::pair<size_t, std::pair<size_t, size_t>>& popTrio)
{
  double f = getParamterValue("f_p2_p3");
}
