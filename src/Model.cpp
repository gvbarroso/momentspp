/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 08/05/2023
 *
 */


#include "Model.hpp"

void Model::fireParameterChanged(const bpp::ParameterList& params)
{
  assert(data_ != nullptr);

  matchParametersValues(params);
  updateEpochs_(params); // updates transitionMatrix_ within each epoch
  computeExpectedSumStats();
  computeCompositeLogLikelihood_(); // e.g. for each rec. bin
}

void Model::computeExpectedSumStats()
{
  expected_ = epochs_[0]->getSteadyState(); // resets moments to the "deep past"

  for(size_t i = 1; i < epochs_.size() - 1; ++i) // epochs are sorted from past to present
  {
    epochs_[i]->transferStatistics(expected_); // copying values into epoch i according to population ancestry
    epochs_[i]->computeExpectedSumStats(expected_); // trickling moments down epochs
  }

  // final epoch
  epochs_.back()->transferStatistics(expected_);
  epochs_.back()->computeExpectedSumStats(expected_);
  epochs_.back()->updateMoments(expected_); // updates inside sslib
}

void Model::printAliasedMoments(std::ostream& stream)
{
  std::vector<std::shared_ptr<Moment>> tmp = epochs_.back()->getSslib().getBasis();

  for(auto& m : tmp)
    stream << m->getName() << " = " << m->getValue() << "\n";
}

void Model::updateEpochs_(const bpp::ParameterList& params)
{
  for(auto it = std::begin(epochs_); it != std::end(epochs_); ++it)
    (*it)->fireParameterChanged(params);
}

void Model::computeCompositeLogLikelihood_()
{
  double cll = 0.;
  /*for(auto it = std::begin(recBins_); it != std::end(recBins_); ++it)
  {
    Eigen::VectorXd obsMeans = data_->getY();
    Eigen::MatrixXd obsCovarMat = data_->getCovarMatrix();
    cll += det(2*covarMat)^(-1/2) * exp(-1/2 * (expected_ - means).transpose() * covarMat^(-1)*(expected_ - means);
  }*/

  compLogLikelihood_ = cll;
}

// this method defines the relationships among moments from different epochs (w.r.t population indices)
// when pop has 2 ancestors, pick one of them, then apply Admixture as if it were a pulse
void Model::linkMoments_()
{
  for(size_t i = 1; i < epochs_.size(); ++i) // for each epoch starting from the 2nd
  {
    std::cout << "epoch " << i << std::endl;
    // for each moment in focal epoch, set "parent" in previous epoch using population ancestry
    for(auto it = std::begin(epochs_[i]->getMoments()); it != std::end(epochs_[i]->getMoments()); ++it)
    {
      if((*it)->getPrefix() == "DD")
      {
        // indices of populations in parental DD_** moment
        size_t prevP1 = epochs_[i - 1]->getSslib().getNumPops(); // inits to out-of-bounds
        size_t prevP2 = epochs_[i - 1]->getSslib().getNumPops(); // inits to out-of-bounds

        size_t focalP1 = (*it)->getPopIndices()[0];
        size_t p1LeftParentId = epochs_[i]->fetchPop(focalP1)->getLeftParent()->getId();
        size_t p1RightParentId = epochs_[i]->fetchPop(focalP1)->getRightParent()->getId();

        if(p1LeftParentId == p1RightParentId) // population [carry-forward / split] between epochs
          prevP1 = p1LeftParentId;

        else
          prevP1 = p1RightParentId;

        size_t focalP2 = (*it)->getPopIndices()[1];
        size_t p2LeftParentId = epochs_[i]->fetchPop(focalP2)->getLeftParent()->getId();
        size_t p2RightParentId = epochs_[i]->fetchPop(focalP2)->getRightParent()->getId();

        if(p2LeftParentId == p2RightParentId) // population [carry-forward / split] between epochs
          prevP2 = p2LeftParentId;

        else
          prevP2 = p2RightParentId;

        size_t idx = epochs_[i - 1]->getSslib().findCompressedIndex(epochs_[i - 1]->getSslib().findDdIndex(prevP1, prevP2));
        (*it)->setParent(epochs_[i - 1]->getSslib().getBasis()[idx]);
      }

      else if((*it)->getPrefix() == "Dz")
      {
        // indices of populations in parental Dz*** moment
        size_t prevP1 = epochs_[i - 1]->getSslib().getNumPops(); // inits to out-of-bounds
        size_t prevP2 = epochs_[i - 1]->getSslib().getNumPops(); // inits to out-of-bounds
        size_t prevP3 = epochs_[i - 1]->getSslib().getNumPops(); // inits to out-of-bounds

        size_t focalP1 = (*it)->getPopIndices()[0];
        size_t p1LeftParentId = epochs_[i]->fetchPop(focalP1)->getLeftParent()->getId();
        size_t p1RightParentId = epochs_[i]->fetchPop(focalP1)->getRightParent()->getId();

        if(p1LeftParentId == p1RightParentId) // population [carry-forward / split] between epochs
          prevP1 = p1LeftParentId;

        else
          prevP1 = p1RightParentId;

        size_t focalP2 = (*it)->getPopIndices()[1];
        size_t p2LeftParentId = epochs_[i]->fetchPop(focalP2)->getLeftParent()->getId();
        size_t p2RightParentId = epochs_[i]->fetchPop(focalP2)->getRightParent()->getId();

        if(p2LeftParentId == p2RightParentId) // population [carry-forward / split] between epochs
          prevP2 = p2LeftParentId;

        else
          prevP2 = p2RightParentId;

        size_t focalP3 = (*it)->getPopIndices()[2];
        size_t p3LeftParentId = epochs_[i]->fetchPop(focalP3)->getLeftParent()->getId();
        size_t p3RightParentId = epochs_[i]->fetchPop(focalP3)->getRightParent()->getId();

        if(p3LeftParentId == p3RightParentId) // population [carry-forward / split] between epochs
          prevP3 = p3LeftParentId;

        else
          prevP3 = p3RightParentId;

        size_t idx = epochs_[i - 1]->getSslib().findCompressedIndex(epochs_[i - 1]->getSslib().findDzIndex(prevP1, prevP2, prevP3));
        (*it)->setParent(epochs_[i - 1]->getSslib().getBasis()[idx]);
      }

      else if((*it)->getPrefix() == "H")
      {
        // indices of populations in parental H_** moment
        size_t prevP1 = epochs_[i - 1]->getSslib().getNumPops(); // inits to out-of-bounds
        size_t prevP2 = epochs_[i - 1]->getSslib().getNumPops(); // inits to out-of-bounds

        size_t focalP1 = (*it)->getPopIndices()[0];
        size_t p1LeftParentId = epochs_[i]->fetchPop(focalP1)->getLeftParent()->getId();
        size_t p1RightParentId = epochs_[i]->fetchPop(focalP1)->getRightParent()->getId();

        if(p1LeftParentId == p1RightParentId) // population [carry-forward / split] between epochs
          prevP1 = p1LeftParentId;

        else // focalP1 is formed by admixture of two populations in previous epoch
          prevP1 = p1RightParentId;

        size_t focalP2 = (*it)->getPopIndices()[1];
        size_t p2LeftParentId = epochs_[i]->fetchPop(focalP2)->getLeftParent()->getId();
        size_t p2RightParentId = epochs_[i]->fetchPop(focalP2)->getRightParent()->getId();

        if(p2LeftParentId == p2RightParentId) // population [carry-forward / split] between epochs
          prevP2 = p2LeftParentId;

        else
          prevP2 = p2RightParentId;

        size_t idx = epochs_[i - 1]->getSslib().findCompressedIndex(epochs_[i - 1]->getSslib().findHetIndex(prevP1, prevP2));
        (*it)->setParent(epochs_[i - 1]->getSslib().getBasis()[idx]);
      }

      else if((*it)->getPrefix() == "pi2")
      {
        // indices of populations in parental pi2**** moment
        size_t prevP1 = epochs_[i - 1]->getSslib().getNumPops(); // inits to out-of-bounds
        size_t prevP2 = epochs_[i - 1]->getSslib().getNumPops(); // inits to out-of-bounds
        size_t prevP3 = epochs_[i - 1]->getSslib().getNumPops(); // inits to out-of-bounds
        size_t prevP4 = epochs_[i - 1]->getSslib().getNumPops(); // inits to out-of-bounds

        size_t focalP1 = (*it)->getPopIndices()[0];
        size_t p1LeftParentId = epochs_[i]->fetchPop(focalP1)->getLeftParent()->getId();
        size_t p1RightParentId = epochs_[i]->fetchPop(focalP1)->getRightParent()->getId();

        if(p1LeftParentId == p1RightParentId) // population [carry-forward / split] between epochs
          prevP1 = p1LeftParentId;

        else
          prevP1 = p1RightParentId;

        size_t focalP2 = (*it)->getPopIndices()[1];
        size_t p2LeftParentId = epochs_[i]->fetchPop(focalP2)->getLeftParent()->getId();
        size_t p2RightParentId = epochs_[i]->fetchPop(focalP2)->getRightParent()->getId();

        if(p2LeftParentId == p2RightParentId) // population [carry-forward / split] between epochs
          prevP2 = p2LeftParentId;

        else
          prevP2 = p2RightParentId;

        size_t focalP3 = (*it)->getPopIndices()[2];
        size_t p3LeftParentId = epochs_[i]->fetchPop(focalP3)->getLeftParent()->getId();
        size_t p3RightParentId = epochs_[i]->fetchPop(focalP3)->getRightParent()->getId();

        if(p3LeftParentId == p3RightParentId) // population [carry-forward / split] between epochs
          prevP3 = p3LeftParentId;

        else
          prevP3 = p3RightParentId;

        size_t focalP4 = (*it)->getPopIndices()[3];
        size_t p4LeftParentId = epochs_[i]->fetchPop(focalP4)->getLeftParent()->getId();
        size_t p4RightParentId = epochs_[i]->fetchPop(focalP4)->getRightParent()->getId();

        if(p4LeftParentId == p4RightParentId) // population [carry-forward / split] between epochs
          prevP4 = p4LeftParentId;

        else
          prevP4 = p4RightParentId;

        size_t idx = epochs_[i - 1]->getSslib().findCompressedIndex(epochs_[i - 1]->getSslib().findPi2Index(prevP1, prevP2, prevP3, prevP4));
        (*it)->setParent(epochs_[i - 1]->getSslib().getBasis()[idx]);
      }

      else if((*it)->getPrefix() == "I")
        (*it)->setParent(epochs_[i - 1]->getSslib().getDummyMomentCompressed());

      (*it)->printAttributes(std::cout);
    } // ends loop over moments
  } // ends loop over epochs
}
