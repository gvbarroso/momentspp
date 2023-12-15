/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 12/10/2023
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
    epochs_[i]->transferStatistics(expected_); // copying values from epoch i-1 into epoch i according to population ancestry
    epochs_[i]->computeExpectedSumStats(expected_); // trickling moments down epochs
  }

  if(epochs_.size() > 1) // final epoch
  {
    epochs_.back()->transferStatistics(expected_);
    epochs_.back()->computeExpectedSumStats(expected_);
    epochs_.back()->updateMoments(expected_); // updates inside sslib
  }
}

void Model::printAliasedMoments(std::ostream& stream)
{
  // prints expectations for the last (most recent) epoch
  std::vector<std::shared_ptr<Moment>> tmp = epochs_.back()->getSslib().getBasis();

  for(auto& m : tmp)
    stream << std::setprecision(20) << m->getName() << " = " << m->getValue() << "\n";
}

void Model::updateEpochs_(const bpp::ParameterList& params)
{
  for(auto it = std::begin(epochs_); it != std::end(epochs_); ++it)
    (*it)->fireParameterChanged(params);
}

void Model::computeCompositeLogLikelihood_()
{
  assert(data_ != nullptr);

  double cll = 0.;
  /*for(auto it = std::begin(recBins_); it != std::end(recBins_); ++it)
  {
    Eigen::VectorXd obsMeans = data_->getY();
    Eigen::MatrixXd obsCovarMat = data_->getCovarMatrix();
    cll += det(2*covarMat)^(-1/2) * exp(-1/2 * (expected_ - means).transpose() * covarMat^(-1)*(expected_ - means);
  }*/

  compLogLikelihood_ = cll;
}

void Model::compressParameters(bool aliasOverEpochs, bool aliasOverPops)
{
  // epoch[0] should never be 1-generation only (ie, an "Admixture epoch")
  // hence it should always have a full set of parameters, ie, including 'u_*', 'r_*', and 's_*'

  // the order of the following two aliasing dimensions (over populations and over epochs) matters!
  // as implemented, we should first go over populations

  // alias u, r and s among populations from the same epoch IFF they have identical (starting) values
  if(aliasOverPops)
  {
    std::cout << "Aliasing parameters over populations.\n";

    for(size_t i = 0; i < epochs_.size(); ++i)
    {
      for(size_t j = 0; j < epochs_[i]->getNumPops(); ++j)
      {
        size_t jd = epochs_[i]->getPops()[j]->getId();
        std::string rj = "r_" + bpp::TextTools::toString(jd);
        std::string uj = "u_" + bpp::TextTools::toString(jd);
        std::string sj = "s_" + bpp::TextTools::toString(jd);

        for(size_t k = j + 1; k < epochs_[i]->getNumPops(); ++k)
        {
          size_t kd = epochs_[i]->getPops()[k]->getId();
          std::string rk = "r_" + bpp::TextTools::toString(kd);
          std::string uk = "u_" + bpp::TextTools::toString(kd);
          std::string sk = "s_" + bpp::TextTools::toString(kd);

          if(epochs_[i]->getParameter(rj).getValue() == epochs_[i]->getParameter(rk).getValue())
          {
            epochs_[i]->aliasParameters(rj, rk);
            aliasParameters(epochs_[i]->getName() + "." + rj, epochs_[i]->getName() + "." + rk);
          }

          if(epochs_[i]->getParameter(uj).getValue() == epochs_[i]->getParameter(uk).getValue())
          {
            epochs_[i]->aliasParameters(uj, uk);
            aliasParameters(epochs_[i]->getName() + "." + uj, epochs_[i]->getName() + "." + uk);
          }

          if(epochs_[i]->getParameter(sj).getValue() == epochs_[i]->getParameter(sk).getValue())
          {
            epochs_[i]->aliasParameters(sj, sk);
            aliasParameters(epochs_[i]->getName() + "." + sj, epochs_[i]->getName() + "." + sk);
          }
        }
      }
    }
  }

  if(aliasOverEpochs)
  {
    std::cout << "Aliasing parameters over epochs.\n";

    for(size_t i = 1; i < epochs_.size(); ++i)
    {
      for(size_t j = 0; j < epochs_[i]->getNumPops(); ++j)
      {
        size_t jd = epochs_[i]->getPops()[j]->getId();
        size_t pd = epochs_[i]->getPops()[j]->getLeftParent()->getId(); // WARNING getLeftParent() [dangerous when there is admixture]

        // rates from population jd
        std::string rj = "r_" + bpp::TextTools::toString(jd);
        std::string uj = "u_" + bpp::TextTools::toString(jd);
        std::string sj = "s_" + bpp::TextTools::toString(jd);

        // rates based on its parental pop id
        std::string rp = "r_" + bpp::TextTools::toString(pd);
        std::string up = "u_" + bpp::TextTools::toString(pd);
        std::string sp = "s_" + bpp::TextTools::toString(pd);

        // only alias if populations share NAME
        if(epochs_[i]->hasIndependentParameter(rj) && epochs_[i]->getPops()[j]->getName() == epochs_[i]->getPops()[j]->getLeftParent()->getName())
          aliasParameters(epochs_[i - 1]->getName() + "." + rp, epochs_[i]->getName() + "." + rj);

        if(epochs_[i]->hasIndependentParameter(uj) && epochs_[i]->getPops()[j]->getName() == epochs_[i]->getPops()[j]->getLeftParent()->getName())
          aliasParameters(epochs_[i - 1]->getName() + "." + up, epochs_[i]->getName() + "." + uj);

        if(epochs_[i]->hasIndependentParameter(sj) && epochs_[i]->getPops()[j]->getName() == epochs_[i]->getPops()[j]->getLeftParent()->getName())
          aliasParameters(epochs_[i - 1]->getName() + "." + sp, epochs_[i]->getName() + "." + sj);
      }
    }
  }
}

// this method defines the relationships among moments from different epochs (w.r.t population indices)
// when pop has 2 ancestors, pick one of them, then apply Admixture as if it were a pulse
void Model::linkMoments_()
{
  for(size_t i = 1; i < epochs_.size(); ++i) // for each epoch starting from the 2nd
  {
    // for each moment in focal epoch, set "parent" in previous epoch using population ancestry
    for(auto it = std::begin(epochs_[i]->getBasis()); it != std::end(epochs_[i]->getBasis()); ++it)
    {
      std::vector<size_t> factorIds = (*it)->getFactorIndices();

      for(size_t j = 0; j < factorIds.size(); ++j)
      {
        size_t prevFactorPop = epochs_[i - 1]->getSslib().getNumPops(); // inits to out-of-bounds
        size_t focalFactorPop = factorIds[j];
        size_t factorPopLeftParentId = epochs_[i]->fetchPop(focalFactorPop)->getLeftParent()->getId();
        size_t factorPopRightParentId = epochs_[i]->fetchPop(focalFactorPop)->getRightParent()->getId();

        if(factorPopLeftParentId == factorPopRightParentId) // population [carry-forward / split] between epochs
          prevFactorPop = factorPopLeftParentId;

        else
          prevFactorPop = factorPopRightParentId;

        factorIds[j] = prevFactorPop; // replaces
      }

      if((*it)->getPrefix() == "DD")
      {
        // indices of populations in parental DD_** moment, inits to out-of-bounds
        size_t prevP1 = epochs_[i - 1]->getSslib().getNumPops();
        size_t prevP2 = epochs_[i - 1]->getSslib().getNumPops();

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

        std::vector<size_t> popIds = { prevP1, prevP2 };
        size_t idx = epochs_[i - 1]->getSslib().findCompressedIndex(epochs_[i - 1]->getSslib().getMoment("DD", popIds, factorIds));
        (*it)->setParent(epochs_[i - 1]->getSslib().getBasis()[idx]);
      }

      #ifdef NAKED_D
      else if((*it)->getPrefix() == "D")
      {
        // indices of populations in parental D_* moment, inits to out-of-bounds
        size_t prevP = epochs_[i - 1]->getSslib().getNumPops();

        size_t focalP = (*it)->getPopIndices()[0];
        size_t pLeftParentId = epochs_[i]->fetchPop(focalP)->getLeftParent()->getId();
        size_t pRightParentId = epochs_[i]->fetchPop(focalP)->getRightParent()->getId();

        if(pLeftParentId == pRightParentId) // population [carry-forward / split] between epochs
          prevP = pLeftParentId;

        else
          prevP = pRightParentId;

        std::vector<size_t>  popIds = { prevP };
        size_t idx = epochs_[i - 1]->getSslib().findCompressedIndex(epochs_[i - 1]->getSslib().getMoment("D", popIds, factorIds));
        (*it)->setParent(epochs_[i - 1]->getSslib().getBasis()[idx]);
      }
      #endif

      else if((*it)->getPrefix() == "Dr")
      {
        // indices of populations in parental Dr*** moment, inits to out-of-bounds
        size_t prevP1 = epochs_[i - 1]->getSslib().getNumPops();
        size_t prevP2 = epochs_[i - 1]->getSslib().getNumPops();

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

        std::vector<size_t> popIds = { prevP1, prevP2 };
        size_t idx = epochs_[i - 1]->getSslib().findCompressedIndex(epochs_[i - 1]->getSslib().getMoment("Dr", popIds, factorIds));
        (*it)->setParent(epochs_[i - 1]->getSslib().getBasis()[idx]);
      }

      else if((*it)->getPrefix() == "Hl")
      {
        // indices of populations in parental Hl_** moment
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

        std::vector<size_t> popIds = { prevP1, prevP2 };
        size_t idx = epochs_[i - 1]->getSslib().findCompressedIndex(epochs_[i - 1]->getSslib().getMoment("Hl", popIds, factorIds));
        (*it)->setParent(epochs_[i - 1]->getSslib().getBasis()[idx]);
      }

      else if((*it)->getPrefix() == "Hr")
      {
        // indices of populations in parental Hr_** moment
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

        std::vector<size_t> popIds = { prevP1, prevP2 };
        size_t idx = epochs_[i - 1]->getSslib().findCompressedIndex(epochs_[i - 1]->getSslib().getMoment("Hr", popIds, factorIds));
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

        std::vector<size_t> popIds = { prevP1, prevP2, prevP3, prevP4 };
        size_t idx = epochs_[i - 1]->getSslib().findCompressedIndex(epochs_[i - 1]->getSslib().getMoment("pi2", popIds, factorIds));
        (*it)->setParent(epochs_[i - 1]->getSslib().getBasis()[idx]);
      }

      else if((*it)->getPrefix() == "I")
        (*it)->setParent(epochs_[i - 1]->getSslib().getMoment("I"));

    } // ends loop over moments
  } // ends loop over epochs
}
