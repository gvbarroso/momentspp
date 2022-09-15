/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 15/09/2022
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

  // we use pops map (of splits and admixtures) to change vector of sum stats between each epoch (copy from i to i + 1)
  for(size_t i = 0; i < epochs_.size() - 1; ++i) // epochs must be sorted from past to present
  {
    epochs_[i]->computeExpectedSumStats(expected_); // trickling sum stats vector down the epochs
    epochs_[i]->copyStatsToMap(expected_); // NOTE can be replaced later by simply finding the indices and copying stats directly into Y

    // WARNING expected_ potentially changes size between epochs

    // for each population in next epoch i + 1
    for(auto it = std::begin(epochs_[i + 1]->getPopsMap()); it != std::end(epochs_[i + 1]->getPopsMap()); ++it)
    {
      size_t idx = (*it).first; // pop index in next epoch
      size_t p1 = (*it).second.first; // parent pop in current epoch
      size_t p2 = (*it).second.second; // parent pop in current epoch

      if(p1 == p2) // parents are the same
      {
        if(idx == p1) // stays the same pop between epochs
        {
          epochs_[i]->getStatsMap();
        }

        else // population split, also dealt with by copying sum stats
        {

        }
      }

      else // parents are different, admixture scenario
        popAdmix_();
    }
  }

  epochs_.back()->computeExpectedSumStats(expected_); // final epoch
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
