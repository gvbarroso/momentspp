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

  for(size_t i = 0; i < epochs_.size() - 1; ++i) // epochs must be sorted from past to present
  {
    epochs_[i]->computeExpectedSumStats(expected_); // trickling sum stats vector down the epochs
    epochs_[i]->copyStatsToMap(expected_); // helper for bookkeeping, might do something better later

    // we use pops map (with info w.r.t splits and admixtures) to change vector of sum stats between each epoch (copy from i to i + 1 following pop indices)
    for(auto itP = std::begin(epochs_[i + 1]->getPopsMap()); itP != std::end(epochs_[i + 1]->getPopsMap()); ++itP) // for each population in next epoch i + 1
    {
      auto tmp = epochs_[i]->getStatsMap(); // lazy copy, do first, do better later

      size_t idx = itP->first; // pop index in next epoch
      size_t p1 = itP->second.first; // parent pop in current epoch
      size_t p2 = itP->second.second; // parent pop in current epoch

      if(p1 == p2)
      {
        if(idx == p1) // pop simply carries over to the next epoch
        {
          // searches for stats contaning 'p1' index in epoch_[i], copy their values to stats containing idx in epoch_[i + 1]
          for(auto itS = std::begin(epochs_[i]->getStatsMap()); itS != std::end(epochs_[i]->getStatsMap()); ++itS)
          {
            std::string mom = (*itS).first;
            if(epochs_[i]->getSslib().hasPopIndex(mom, bpp::TextTools::toString(p1))) // if moment has p1 in any pop index
              epochs_[i + 1]->getSslib().setMomValue(mom, (*itS).second); // moment names are identical between epochs in case idx == p1
          }
        }

        else // happens when there is a pop split, eg, (1,(1,1)) ; (2,(1,1)) , meaning p1 == p2 will be visited twice
        {
          // searches for stats contaning 'p1' index in epoch_[i], copy their values to stats containing idx in epoch_[i + 1]
          for(auto itS = std::begin(epochs_[i]->getStatsMap()); itS != std::end(epochs_[i]->getStatsMap()); ++itS)
          {
            std::string mom = (*itS).first;
            if(epochs_[i]->getSslib().hasPopIndex(mom, bpp::TextTools::toString(p1))) // if moment has p1 in any pop index
            {
              std::string momNext = ; // TODO build moment with same name except idx instead of p1
              epochs_[i + 1]->getSslib().setMomValue(mom, (*itS).second); // moment names are identical between epochs in case idx == p1
            }
          }
        }
      }

      else // parents are different, admixture scenario
        popAdmix_(idx, p1, p2, i); // i (or something else) so we can track the epoch down to popAdmix_() method
    }

    expected_ = epochs_[i + 1]->fetchYvec();
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
