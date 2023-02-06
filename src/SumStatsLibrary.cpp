/*
 * Authors: Gustavo V. Barroso
 * Created: 05/08/2022
 * Last modified: 06/02/2023
 *
 */

#include <ios>

#include "SumStatsLibrary.hpp"

Eigen::VectorXd SumStatsLibrary::fetchYvec()
{
  Eigen::VectorXd y(moments_.size());

  for(size_t i = 0; i < moments_.size(); ++i)
    y(i) = moments_[i]->getValue();

  return y;
}

void SumStatsLibrary::printMoments(std::ostream& stream)
{
  for(size_t i = 0; i < moments_.size(); ++i)
    moments_[i]->printAttributes(stream);
}

void SumStatsLibrary::initMoments_()
{
  moments_.reserve(getNumStats());

  // first pass to include DD, Dz, H statistics
  for(auto itI = std::begin(popIndices_); itI != std::end(popIndices_); ++itI)
  {
    for(auto itJ = std::begin(popIndices_); itJ != std::end(popIndices_); ++itJ)
    {
      moments_.emplace_back(std::make_shared<DdMoment>("DD_" + asString(*itI) + "_" + asString(*itJ) + "_X", 0.));

      moments_.emplace_back(std::make_shared<HetMoment>("H_" + asString(*itI) + "_" + asString(*itJ) + "_A", 0., true, false)); // H_ii p(1-p)
      moments_.emplace_back(std::make_shared<HetMoment>("H_" + asString(*itI) + "_" + asString(*itJ) + "_B", 0., false, false)); // H_ii (1-p)p

      for(auto itK = std::begin(popIndices_); itK != std::end(popIndices_); ++itK)
        moments_.emplace_back(std::make_shared<DzMoment>("Dz_" + asString(*itI) + "_" + asString(*itJ) + "_" + asString(*itK) + "_X", 0.));
    }
  }

  // second pass to include pi2 statistics, which need pointers to H statistics
  for(auto itI = std::begin(popIndices_); itI != std::end(popIndices_); ++itI)
  {
    for(auto itJ = std::begin(popIndices_); itJ != std::end(popIndices_); ++itJ)
    {
      for(auto itK = std::begin(popIndices_); itK != std::end(popIndices_); ++itK)
      {
        for(auto itL = std::begin(popIndices_); itL != std::end(popIndices_); ++itL)
        {
          std::shared_ptr<HetMoment> leftA;
          std::shared_ptr<HetMoment> leftB;

          std::shared_ptr<HetMoment> rightA;
          std::shared_ptr<HetMoment> rightB;

          // naive search for corresponding H stats because moments_ has not been finalized yet (indices should not be used)
          for(auto itMom = std::begin(moments_); itMom != std::end(moments_); ++itMom)
          {
            // we can't use else if() in the conditionals below because it leads to missed hits when K and L are the same as I and J indices
            if((*itMom)->getName() == "H_" + asString(*itI) + "_" + asString(*itJ) + "_A") // left locus, p(1-p)
            {
              HetMoment* tmp = dynamic_cast<HetMoment*>(itMom->get());
              leftA.reset(tmp);
            }

            if((*itMom)->getName() == "H_" + asString(*itI) + "_" + asString(*itJ) + "_B") // left locus, (1-p)p
            {
              HetMoment* tmp = dynamic_cast<HetMoment*>(itMom->get());
              leftB.reset(tmp);
            }

            if((*itMom)->getName() == "H_" + asString(*itK) + "_" + asString(*itL) + "_A") // right locus, p(1-p)
            {
              HetMoment* tmp = dynamic_cast<HetMoment*>(itMom->get());
              rightA.reset(tmp);
            }

            if((*itMom)->getName() == "H_" + asString(*itK) + "_" + asString(*itL) + "_B") // right locus, (1-p)p
            {
              HetMoment* tmp = dynamic_cast<HetMoment*>(itMom->get());
              rightB.reset(tmp);
            }
          }

          moments_.emplace_back(std::make_shared<Pi2Moment>("pi2_" + asString(*itI) + "_" + asString(*itJ) + "_" + asString(*itK) + "_" + asString(*itL) + "_A", 0., leftA, rightA));
          moments_.emplace_back(std::make_shared<Pi2Moment>("pi2_" + asString(*itI) + "_" + asString(*itJ) + "_" + asString(*itK) + "_" + asString(*itL) + "_B", 0., leftA, rightB));
          moments_.emplace_back(std::make_shared<Pi2Moment>("pi2_" + asString(*itI) + "_" + asString(*itJ) + "_" + asString(*itK) + "_" + asString(*itL) + "_C", 0., leftB, rightA));
          moments_.emplace_back(std::make_shared<Pi2Moment>("pi2_" + asString(*itI) + "_" + asString(*itJ) + "_" + asString(*itK) + "_" + asString(*itL) + "_D", 0., leftB, rightB));
        }
      }
    }
  }

  // adds Dummy Moment lexicographically after H_ stats to convert into a homogeneous system (see Mutation::setUpMatrices_())
  moments_.emplace_back(std::make_shared<Moment>("I", 1.));

  // this determines the ascending lexicographical order of stats in the rows of transition matrices
  std::sort(std::begin(moments_), std::end(moments_), [](std::shared_ptr<Moment> a, std::shared_ptr<Moment> b) { return a->getName() < b->getName(); } );

  for(size_t i = 0; i < moments_.size(); ++i)
    moments_[i]->setPosition(i);

  printMoments(std::cout);
}

// proposal: argument tells in which Populations (id) selection acts on the derived allele of left locus?
void SumStatsLibrary::compressBasis_(const std::vector<size_t>& selectedPopIds)
{
  // ?
  //for(size_t i = 0; i < selectedPopIds.size(); ++i)
  //{
  //  if()
  //  areHetsPermuted_ = false;
  //  arePi2sPermuted_ = true;
  //}
}
