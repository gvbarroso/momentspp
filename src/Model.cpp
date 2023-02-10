/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 09/02/2023
 *
 */


#include "Model.hpp"

void Model::fireParameterChanged(const bpp::ParameterList& params)
{
  matchParametersValues(params);
  updateEpochs_(params); // updates transitionMatrix_ within each epoch

  computeExpectedSumStats();
  computeCompositeLogLikelihood_(data_.getObsY(), data_.getCovarMatrix()); // e.g. for each rec. binx
}

void Model::updateEpochs_(const bpp::ParameterList& params)
{
  for(auto it = std::begin(epochs_); it != std::end(epochs_); ++it)
    (*it)->fireParameterChanged(params);
}

void Model::computeExpectedSumStats()
{
  std::ofstream stats;

  expected_ = epochs_[0]->getSteadyState(); // resets moments to the "deep past"

  stats.open(name_ + "_steady_state.txt");
  epochs_[0]->getSslib().printMoments(stats);
  stats.close();

  for(size_t i = 0; i < epochs_.size() - 1; ++i) // epochs are sorted from past to present
  {
    epochs_[i]->computeExpectedSumStats(expected_); // trickling moments down epochs (pass expected_ by ref)
    epochs_[i + 1]->transferStatistics(expected_); // copying values into epoch i + 1 according to population ancestry (pass expected_ by ref)
  }

  epochs_.back()->computeExpectedSumStats(expected_); // final epoch (out of the for loop due to "i+1" access there)
  epochs_.back()->updateMoments(expected_);

  stats.open(name_ + "_final_moments.txt");
  epochs_.back()->getSslib().printMoments(stats);
  stats.close();
}

void Model::aliasMoments()
{
  // epochs have their own Population container, and Population objects tell if the left derived allele is suject to selection in that epoch/pop
  for(size_t i = 0; i < epochs_.size(); ++i)
  {
    epochs_[i]->getSslib().aliasMoments(epochs_[i]->fetchSelectedPopIds());
    auto tmp = epochs_[i]->getSslib().fetchCompressedBasis();

    std::cout << epochs_[i]->getName() << "'s Moments w/ unique expectations:\n";
    for(auto& m : tmp)
      std::cout << m->getName() << " = " << m->getValue() * (m->getAliases().size() + 1) << "\n";
  }
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

// this method defines the relationships among moments from different epochs (w.r.t population indices)
void Model::linkMoments_()
{
  for(size_t i = 1; i < epochs_.size(); ++i) // for each epoch starting from the 2nd
  {
    // for each moment in focal epoch, set "parents" in previous epoch using input population ancestry
    for(auto it = std::begin(epochs_[i]->getMoments()); it != std::end(epochs_[i]->getMoments()); ++it)
    {
      if((*it)->getPrefix() == "DD")
      {
        // indices of populations in parental DD_** moment
        size_t prevP1 = epochs_[i - 1]->getSslib().getNumPops(); // inits to out-of-bounds
        size_t prevP2 = epochs_[i - 1]->getSslib().getNumPops(); // inits to out-of-bounds

        size_t focalP1 = (*it)->getPopIndices()[0]; // i in DD_i_*
        size_t p1LeftParentId = epochs_[i]->getPops().at(focalP1)->getLeftParent()->getId();
        size_t p1RightParentId = epochs_[i]->getPops().at(focalP1)->getRightParent()->getId();

        if(p1LeftParentId == p1RightParentId) // population [carry-forward / split] between epochs
          prevP1 = p1LeftParentId;

        else
          throw bpp::Exception("Model::TODO->include admixture cases");

        size_t focalP2 = (*it)->getPopIndices()[1]; // j in DD_*_j
        size_t p2LeftParentId = epochs_[i]->getPops().at(focalP2)->getLeftParent()->getId();
        size_t p2RightParentId = epochs_[i]->getPops().at(focalP2)->getRightParent()->getId();

        if(p2LeftParentId == p2RightParentId) // population [carry-forward / split] between epochs
          prevP2 = p2LeftParentId;

        else
          throw bpp::Exception("Model::TODO->include admixture cases");

        (*it)->setParent(epochs_[i - 1]->getSslib().getDdMoment(prevP1, prevP2));
      }

      else if((*it)->getPrefix() == "Dz")
      {
        // indices of populations in parental Dz*** moment
        size_t prevP1 = epochs_[i - 1]->getSslib().getNumPops(); // inits to out-of-bounds
        size_t prevP2 = epochs_[i - 1]->getSslib().getNumPops(); // inits to out-of-bounds
        size_t prevP3 = epochs_[i - 1]->getSslib().getNumPops(); // inits to out-of-bounds

        size_t focalP1 = (*it)->getPopIndices()[0]; // i in Dz_i_j_k
        size_t p1LeftParentId = epochs_[i]->getPops().at(focalP1)->getLeftParent()->getId();
        size_t p1RightParentId = epochs_[i]->getPops().at(focalP1)->getRightParent()->getId();

        if(p1LeftParentId == p1RightParentId) // population [carry-forward / split] between epochs
          prevP1 = p1LeftParentId;

        else
          throw bpp::Exception("Model::TODO->include admixture cases");

        size_t focalP2 = (*it)->getPopIndices()[1]; // j in Dz_i_j_k
        size_t p2LeftParentId = epochs_[i]->getPops().at(focalP2)->getLeftParent()->getId();
        size_t p2RightParentId = epochs_[i]->getPops().at(focalP2)->getRightParent()->getId();

        if(p2LeftParentId == p2RightParentId) // population [carry-forward / split] between epochs
          prevP2 = p2LeftParentId;

        else
          throw bpp::Exception("Model::TODO->include admixture cases");

        size_t focalP3 = (*it)->getPopIndices()[2]; // k in Dz_i_j_k
        size_t p3LeftParentId = epochs_[i]->getPops().at(focalP3)->getLeftParent()->getId();
        size_t p3RightParentId = epochs_[i]->getPops().at(focalP3)->getRightParent()->getId();

        if(p3LeftParentId == p3RightParentId) // population [carry-forward / split] between epochs
          prevP3 = p3LeftParentId;

        else
          throw bpp::Exception("Model::TODO->include admixture cases");

        (*it)->setParent(epochs_[i - 1]->getSslib().getDzMoment(prevP1, prevP2, prevP3));
      }

      else if((*it)->getPrefix() == "H")
      {
        std::string suffix = (*it)->getSuffix();

        // indices of populations in parental H_**_* moment
        size_t prevP1 = epochs_[i - 1]->getSslib().getNumPops(); // inits to out-of-bounds
        size_t prevP2 = epochs_[i - 1]->getSslib().getNumPops(); // inits to out-of-bounds

        size_t focalP1 = (*it)->getPopIndices()[0]; // i in H_i_*
        size_t p1LeftParentId = epochs_[i]->getPops().at(focalP1)->getLeftParent()->getId();
        size_t p1RightParentId = epochs_[i]->getPops().at(focalP1)->getRightParent()->getId();

        if(p1LeftParentId == p1RightParentId) // population [carry-forward / split] between epochs
          prevP1 = p1LeftParentId;

        else
          throw bpp::Exception("Model::TODO->include admixture cases");

        size_t focalP2 = (*it)->getPopIndices()[1]; // j in H_*_j
        size_t p2LeftParentId = epochs_[i]->getPops().at(focalP2)->getLeftParent()->getId();
        size_t p2RightParentId = epochs_[i]->getPops().at(focalP2)->getRightParent()->getId();

        if(p2LeftParentId == p2RightParentId) // population [carry-forward / split] between epochs
          prevP2 = p2LeftParentId;

        else
          throw bpp::Exception("Model::TODO->include admixture cases");

        (*it)->setParent(epochs_[i - 1]->getSslib().getHetMoment(prevP1, prevP2, suffix));
      }

      else if((*it)->getPrefix() == "pi2")
      {
        std::string suffix = (*it)->getSuffix();

        // indices of populations in parental pi2**** moment
        size_t prevP1 = epochs_[i - 1]->getSslib().getNumPops(); // inits to out-of-bounds
        size_t prevP2 = epochs_[i - 1]->getSslib().getNumPops(); // inits to out-of-bounds
        size_t prevP3 = epochs_[i - 1]->getSslib().getNumPops(); // inits to out-of-bounds
        size_t prevP4 = epochs_[i - 1]->getSslib().getNumPops(); // inits to out-of-bounds

        size_t focalP1 = (*it)->getPopIndices()[0]; // i in pi2_i_j_k_l
        size_t p1LeftParentId = epochs_[i]->getPops().at(focalP1)->getLeftParent()->getId();
        size_t p1RightParentId = epochs_[i]->getPops().at(focalP1)->getRightParent()->getId();

        if(p1LeftParentId == p1RightParentId) // population [carry-forward / split] between epochs
          prevP1 = p1LeftParentId;

        else
          throw bpp::Exception("Model::TODO->include admixture cases");

        size_t focalP2 = (*it)->getPopIndices()[1]; // second population id of epochs i+1's H_** moment
        size_t p2LeftParentId = epochs_[i]->getPops().at(focalP2)->getLeftParent()->getId();
        size_t p2RightParentId = epochs_[i]->getPops().at(focalP2)->getRightParent()->getId();

        if(p2LeftParentId == p2RightParentId) // population [carry-forward / split] between epochs
          prevP2 = p2LeftParentId;

        else
          throw bpp::Exception("Model::TODO->include admixture cases");

        size_t focalP3 = (*it)->getPopIndices()[2]; // second population id of epochs i+1's H_** moment
        size_t p3LeftParentId = epochs_[i]->getPops().at(focalP3)->getLeftParent()->getId();
        size_t p3RightParentId = epochs_[i]->getPops().at(focalP3)->getRightParent()->getId();

        if(p3LeftParentId == p3RightParentId) // population [carry-forward / split] between epochs
          prevP3 = p3LeftParentId;

        else
          throw bpp::Exception("Model::TODO->include admixture cases");

        size_t focalP4 = (*it)->getPopIndices()[3]; // second population id of epochs i+1's H_** moment
        size_t p4LeftParentId = epochs_[i]->getPops().at(focalP4)->getLeftParent()->getId();
        size_t p4RightParentId = epochs_[i]->getPops().at(focalP4)->getRightParent()->getId();

        if(p4LeftParentId == p4RightParentId) // population [carry-forward / split] between epochs
          prevP4 = p4LeftParentId;

        else
          throw bpp::Exception("Model::TODO->include admixture cases");

        (*it)->setParent(epochs_[i - 1]->getSslib().getPi2Moment(prevP1, prevP2, prevP3, prevP4, suffix));
      }

      else if((*it)->getPrefix() == "I")
        (*it)->setParent(epochs_[i - 1]->getSslib().getDummyMoment());

    } // ends for loop over moments
  } // ends for loop over epochs
}
