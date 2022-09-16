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
    Eigen::VectorXd tmp = epochs_[i + 1]->fetchYvec(); // expected and tmp may have different sizes

    // we use pops map (which has info w.r.t splits and admixtures) to change vector of sum stats between each epoch (copy from i to i + 1 following pop indices path)
    for(auto itP = std::begin(epochs_[i + 1]->getPopsMap()); itP != std::end(epochs_[i + 1]->getPopsMap()); ++itP) // for each population in next epoch i + 1
    {
      size_t id = itP->first; // pop index in next epoch
      size_t p1 = itP->second.first; // parent pop in current epoch
      size_t p2 = itP->second.second; // parent pop in current epoch

      if(p1 == p2) // pop split / carry-over
      {
        popSplit_(epochs_[i], epochs_[i + 1], p1, id);

        expected_ = tmp; // swap
      }

      else // admixture event
        popAdmix_(idx, p1, p2, i); // i (or something else) so we can track the epoch down to popAdmix_() method
    }

    expected_ = epochs_[i + 1]->fetchYvec();
  }

  epochs_.back()->computeExpectedSumStats(expected_); // final epoch (out of the for loop due to "i+1" access there)
}

void Model::popSplit_(std::shard_ptr<Epoch> epochFrom, std::shard_ptr<Epoch> epochTo, size_t popIdFrom, size_t popIdTo)
{
  std::vector<size_t> indicesFrom(0);
  std::vector<size_t> indicesTo(0);

  // TODO pass vectors as non-const ref to find*StatIndices() methods

  // searches for positions of Het stats contaning popIdFrom pop-index in epoch_[i] and popIdTo pop-index in epoch_[i + 1]
  indicesFrom = epochFrom->getSslib().findHetStatIndices(popIdFrom);
  indicesTo = epochTo->getSslib().findHetStatIndices(popIdTo);

  // copy their values to Het stats containing '
  for(size_t j = 0; j < indicesFrom.size(); ++j)
  {
    double value = expected_(0, indicesFrom[j]);
  }

  // searches for positions of Pi2 stats contaning popIdFrom pop-index in epoch_[i] and popIdTo pop-index in epoch_[i + 1]
  indicesFrom = epochFrom->getSslib().findPi2StatIndices(popIdFrom);
  indicesTo = epochTo->getSslib().findPi2StatIndices(popIdTo);

  // copy their values to Pi2 stats containing '
  for(size_t j = 0; j < indicesFrom.size(); ++j)
  {
    double value = expected_(0, indicesFrom[j]);
  }

  // searches for positions of DD stats contaning popIdFrom pop-index in epoch_[i] and popIdTo pop-index in epoch_[i + 1]
  indicesFrom = epochFrom->getSslib().findDDStatIndices(popIdFrom);
  indicesTo = epochTo->getSslib().findDDStatIndices(popIdTo);

  // copy their values to DD stats containing '
  for(size_t j = 0; j < indicesFrom.size(); ++j)
  {
    double value = expected_(0, indicesFrom[j]);
  }

  // searches for positions of Dz stats contaning popIdFrom pop-index in epoch_[i] and popIdTo pop-index in epoch_[i + 1]
  indicesFrom = epochFrom->getSslib().findDzStatIndices(popIdFrom);
  indicesTo = epochTo->getSslib().findDzStatIndices(popIdTo);

  // copy their values to Dz stats containing '
  for(size_t j = 0; j < indicesFrom.size(); ++j)
  {
    double value = expected_(0, indicesFrom[j]);
  }

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
