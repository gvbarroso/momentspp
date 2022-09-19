/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 16/09/2022
 *
 */


#include "Model.hpp"

void Model::fireParameterChanged(const bpp::ParameterList& params)
{
  matchParametersValues(params);
  updateEpochs_(params); // updates transitionMatrix_ within each epoch

  computeExpectedSumStats_();
  computeCompositeLogLikelihood_(data_.getYvec(), data_.getCovarMatrix()); // for each rec. binx
}

void Model::updateEpochs_(const bpp::ParameterList& params)
{
  for(auto it = std::begin(epochs_); it != std::end(epochs_); ++it)
    (*it)->fireParameterChanged(params);
}

void Model::computeExpectedSumStats_()
{
  expected_ = epochs_[0]->getSteadyState(); // resets stats to the "deep past"

  for(size_t i = 0; i < epochs_.size() - 1; ++i) // epochs must be sorted from past to present
  {
    epochs_[i]->computeExpectedSumStats(expected_); // trickling sum stats vector down the epochs (non-const ref)

    // for each statistic in next epoch i + 1
    for(auto itS = std::begin(epochs_[i + 1]->getStatsMap()); itS != std::end(epochs_[i + 1]->getStatsMap()); ++itS)
    {
      std::string momTo = (*itS)->first;
      std::vector<std::string> splitMom = sslib.splitString(momTo, "_");

      if(splitMom[0] == "DD")
      {

      }

      else if(splitMom[0] == "Dz")
      {

      }

      else if(splitMom[0] == "H")
      {
        size_t p1 = splitMom[0][0]; // first population id of epochs_[i + 1] moment
        size_t p1LeftParentId = epochs_[i + 1]->getPops().at(p1)->getLeftParent()->getId();
        size_t p1RightParentId = epochs_[i + 1]->getPops().at(p1)->getRightParent()->getId();

        if(p1LeftParentId == p1RightParentId) // population carry-forward / split between epochs
          p1 = p1LeftParentId;

        size_t p2 = splitMom[0][1]; // first population id of epochs_[i + 1] moment
        size_t p2LeftParentId = epochs_[i + 1]->getPops().at(p2)->getLeftParent()->getId();
        size_t p2RightParentId = epochs_[i + 1]->getPops().at(p2)->getRightParent()->getId();

        if(p2LeftParentId == p2RightParentId) // population carry-forward / split between epochs
          p2 = p2LeftParentId;

        std::string momFrom = "H_" + bpp::TextTools::toString(p1) + bpp::TextTools::toString(p2);
        double cpyVal = epochs_[i]->getStatsMap().at(momFrom).second;

        epochs_[i + 1]->getSslib().setValue(momTo, cpyVal);
      }

      else if(splitMom[0] == "pi2")
      {

      }

        expected_ = tmp; // swap


      else // admixture event
        popAdmix_(idx, p1, p2, i); // i (or something else) so we can track the epoch down to popAdmix_() method
    }

    expected_ = epochs_[i + 1]->fetchYvec();
  }

  epochs_.back()->computeExpectedSumStats(expected_); // final epoch (out of the for loop due to "i+1" access there)
}

void Model::popSplit_(std::shard_ptr<Epoch> epochFrom, std::shard_ptr<Epoch> epochTo, size_t popIdFrom, size_t popIdTo)
{

}

// admixes p2 and p3 (from second)
void Model::popAdmix_(const std::pair<size_t, std::pair<size_t, size_t>>& popTrio)
{
  double f = getParamterValue("f_p2_p3");
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
