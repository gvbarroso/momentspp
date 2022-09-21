/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 20/09/2022
 *
 */


#include "Model.hpp"

void Model::fireParameterChanged(const bpp::ParameterList& params)
{
  matchParametersValues(params);
  updateEpochs_(params); // updates transitionMatrix_ within each epoch

  computeExpectedSumStats_();
  computeCompositeLogLikelihood_(data_.getYvec(), data_.getCovarMatrix()); // e.g. for each rec. binx
}

void Model::updateEpochs_(const bpp::ParameterList& params)
{
  for(auto it = std::begin(epochs_); it != std::end(epochs_); ++it)
    (*it)->fireParameterChanged(params);
}

void Model::computeExpectedSumStats_()
{
  expected_ = epochs_[0]->getSteadyState(); // resets moments to the "deep past"

  for(size_t i = 0; i < epochs_.size() - 1; ++i) // epochs must be sorted from past to present
  {
    epochs_[i]->computeExpectedSumStats(expected_); // trickling moments vector down the epochs (passing expected_ by non-const ref)
    // TODO set values of Moments from epoch i using expected_

    // we want to copy values of Moments in epoch i into corresponding Moments in epoch i + 1 (following population ancestry mapping)
    for(auto itMom = std::begin(epochs_[i + 1]->getMoments()); itMom != std::end(epochs_[i + 1]->getMoments()); ++itMom)
    {
      if(itMom->getPrefix() == "DD")
      {
        size_t prevP1, prevP2 = 0;

        size_t focalP1 = itMom->getPopIndices()[0]; // first population id of epochs i+1's H_** moment
        size_t p1LeftParentId = epochs_[i + 1]->getPops().at(focalP1)->getLeftParent()->getId();
        size_t p1RightParentId = epochs_[i + 1]->getPops().at(focalP1)->getRightParent()->getId();

        if(p1LeftParentId == p1RightParentId) // population [carry-forward / split] between epochs
          prevP1 = p1LeftParentId;

        size_t focalP2 = itMom->getPopIndices()[1]; // second population id of epochs i+1's H_** moment
        size_t p2LeftParentId = epochs_[i + 1]->getPops().at(focalP2)->getLeftParent()->getId();
        size_t p2RightParentId = epochs_[i + 1]->getPops().at(focalP2)->getRightParent()->getId();

        if(p2LeftParentId == p2RightParentId) // population [carry-forward / split] between epochs
          prevP2 = p2LeftParentId;

        double cpyVal = epochs_[i]->getSslib().getDdMoment(prevP1, prevP2).getValue();
        epochs_[i + 1]->getSslib().setHetMomentValue(focalP1, focalP2, cpyVal);
      }

      else if(itMom->getPrefix() == "Dz")
      {
        size_t prevP1, prevP2, prevP3 = 0;

        size_t focalP1 = itMom->getPopIndices()[0]; // first population id of epochs i+1's H_** moment
        size_t p1LeftParentId = epochs_[i + 1]->getPops().at(focalP1)->getLeftParent()->getId();
        size_t p1RightParentId = epochs_[i + 1]->getPops().at(focalP1)->getRightParent()->getId();

        if(p1LeftParentId == p1RightParentId) // population [carry-forward / split] between epochs
          prevP1 = p1LeftParentId;

        size_t focalP2 = itMom->getPopIndices()[1]; // second population id of epochs i+1's H_** moment
        size_t p2LeftParentId = epochs_[i + 1]->getPops().at(focalP2)->getLeftParent()->getId();
        size_t p2RightParentId = epochs_[i + 1]->getPops().at(focalP2)->getRightParent()->getId();

        if(p2LeftParentId == p2RightParentId) // population [carry-forward / split] between epochs
          prevP2 = p2LeftParentId;

        size_t focalP3 = itMom->getPopIndices()[2]; // second population id of epochs i+1's H_** moment
        size_t p3LeftParentId = epochs_[i + 1]->getPops().at(focalP3)->getLeftParent()->getId();
        size_t p3RightParentId = epochs_[i + 1]->getPops().at(focalP3)->getRightParent()->getId();

        if(p3LeftParentId == p3RightParentId) // population [carry-forward / split] between epochs
          prevP3 = p3LeftParentId;

        double cpyVal = epochs_[i]->getSslib().getDzMoment(prevP1, prevP2, prevP3).getValue();
        epochs_[i + 1]->getSslib().setHetMomentValue(focalP1, focalP2, cpyVal);
      }

      else if(itMom->getPrefix() == "H")
      {
        size_t prevP1, prevP2 = 0;

        size_t focalP1 = itMom->getPopIndices()[0]; // first population id of epochs i+1's H_** moment
        size_t p1LeftParentId = epochs_[i + 1]->getPops().at(focalP1)->getLeftParent()->getId();
        size_t p1RightParentId = epochs_[i + 1]->getPops().at(focalP1)->getRightParent()->getId();

        if(p1LeftParentId == p1RightParentId) // population [carry-forward / split] between epochs
          prevP1 = p1LeftParentId;

        size_t focalP2 = itMom->getPopIndices()[1]; // second population id of epochs i+1's H_** moment
        size_t p2LeftParentId = epochs_[i + 1]->getPops().at(focalP2)->getLeftParent()->getId();
        size_t p2RightParentId = epochs_[i + 1]->getPops().at(focalP2)->getRightParent()->getId();

        if(p2LeftParentId == p2RightParentId) // population [carry-forward / split] between epochs
          prevP2 = p2LeftParentId;

        double cpyVal = epochs_[i]->getSslib().getHetMoment(prevP1, prevP2).getValue();
        epochs_[i + 1]->getSslib().setHetMomentValue(focalP1, focalP2, cpyVal);
      }

      else if(itMom->getPrefix() == "pi2")
      {
        size_t prevP1, prevP2, prevP3, prevP4 = 0;

        size_t focalP1 = itMom->getPopIndices()[0]; // first population id of epochs i+1's H_** moment
        size_t p1LeftParentId = epochs_[i + 1]->getPops().at(focalP1)->getLeftParent()->getId();
        size_t p1RightParentId = epochs_[i + 1]->getPops().at(focalP1)->getRightParent()->getId();

        if(p1LeftParentId == p1RightParentId) // population [carry-forward / split] between epochs
          prevP1 = p1LeftParentId;

        size_t focalP2 = itMom->getPopIndices()[1]; // second population id of epochs i+1's H_** moment
        size_t p2LeftParentId = epochs_[i + 1]->getPops().at(focalP2)->getLeftParent()->getId();
        size_t p2RightParentId = epochs_[i + 1]->getPops().at(focalP2)->getRightParent()->getId();

        if(p2LeftParentId == p2RightParentId) // population [carry-forward / split] between epochs
          prevP2 = p2LeftParentId;

        size_t focalP3 = itMom->getPopIndices()[2]; // second population id of epochs i+1's H_** moment
        size_t p3LeftParentId = epochs_[i + 1]->getPops().at(focalP3)->getLeftParent()->getId();
        size_t p3RightParentId = epochs_[i + 1]->getPops().at(focalP3)->getRightParent()->getId();

        if(p3LeftParentId == p3RightParentId) // population [carry-forward / split] between epochs
          prevP3 = p3LeftParentId;

        size_t focalP4 = itMom->getPopIndices()[3]; // second population id of epochs i+1's H_** moment
        size_t p4LeftParentId = epochs_[i + 1]->getPops().at(focalP4)->getLeftParent()->getId();
        size_t p4RightParentId = epochs_[i + 1]->getPops().at(focalP4)->getRightParent()->getId();

        if(p4LeftParentId == p4RightParentId) // population [carry-forward / split] between epochs
          prevP4 = p4LeftParentId;

        double cpyVal = epochs_[i]->getSslib().getPi2Moment(prevP1, prevP2, prevP3, prevP4).getValue();
        epochs_[i + 1]->getSslib().setPi2MomentValue(focalP1, focalP2, focalP3, focalP4, cpyVal);
      }

      expected_ = epochs_[i + 1]->getSslib().fetchYvec(); // swap
    }
  }

  epochs_.back()->computeExpectedSumStats(expected_); // final epoch (out of the for loop due to "i+1" access there)
}

void Model::popAdmix_()
{
  // TODO
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
