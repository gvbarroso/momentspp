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
      {
        moments_.emplace_back(std::make_shared<DzMoment>("Dz_" + asString(*itI) + "_" + asString(*itJ) + "_" + asString(*itK) + "_X", 0.));

        for(auto itL = std::begin(popIndices_); itL != std::end(popIndices_); ++itL)
        {
          moments_.emplace_back(std::make_shared<Pi2Moment>("pi2_" + asString(*itI) + "_" + asString(*itJ) + "_" + asString(*itK) + "_" + asString(*itL) + "_A", 0., nullptr, nullptr));
          moments_.emplace_back(std::make_shared<Pi2Moment>("pi2_" + asString(*itI) + "_" + asString(*itJ) + "_" + asString(*itK) + "_" + asString(*itL) + "_B", 0., nullptr, nullptr));
          moments_.emplace_back(std::make_shared<Pi2Moment>("pi2_" + asString(*itI) + "_" + asString(*itJ) + "_" + asString(*itK) + "_" + asString(*itL) + "_C", 0., nullptr, nullptr));
          moments_.emplace_back(std::make_shared<Pi2Moment>("pi2_" + asString(*itI) + "_" + asString(*itJ) + "_" + asString(*itK) + "_" + asString(*itL) + "_D", 0., nullptr, nullptr));
        }
      }
    }
  }

  // adds Dummy Moment lexicographically after H_ stats to convert into a homogeneous system (see Mutation::setUpMatrices_())
  moments_.emplace_back(std::make_shared<Moment>("I", 1.));

  // determines the ascending lexicographical order of stats in the rows/cols of matrices inside AbstractOperators
  std::sort(std::begin(moments_), std::end(moments_), [=](std::shared_ptr<Moment> a, std::shared_ptr<Moment> b) { return a->getName() < b->getName(); } );

  for(size_t i = 0; i < moments_.size(); ++i)
    moments_[i]->setPosition(i);

  linkPi2HetStats_();
  printMoments(std::cout);
}

// for each Pi2Moment, sets the two pointers corresponding to HetMoments (left and right loci)
void SumStatsLibrary::linkPi2HetStats_()
{
  for(auto itMom = std::begin(moments_); itMom != std::end(moments_); ++itMom)
  {
    auto tmpPi2 = std::dynamic_pointer_cast<Pi2Moment>(*itMom);

    if(tmpPi2 != nullptr)
    {
      size_t p1 = tmpPi2->getPopIndices()[0];
      size_t p2 = tmpPi2->getPopIndices()[1];
      size_t p3 = tmpPi2->getPopIndices()[2];
      size_t p4 = tmpPi2->getPopIndices()[3];

      std::string suffix = tmpPi2->getSuffix();

      if(suffix == "A")
      {
        auto tmpHetLeft = std::dynamic_pointer_cast<HetMoment>(getHetMoment(p1, p2, "A"));
        auto tmpHetRight = std::dynamic_pointer_cast<HetMoment>(getHetMoment(p3, p4, "A"));
        tmpPi2->setLeftHetStat(tmpHetLeft);
        tmpPi2->setRightHetStat(tmpHetRight);
      }

      else if(suffix == "B")
      {
        auto tmpHetLeft = std::dynamic_pointer_cast<HetMoment>(getHetMoment(p1, p2, "A"));
        auto tmpHetRight = std::dynamic_pointer_cast<HetMoment>(getHetMoment(p3, p4, "B"));
        tmpPi2->setLeftHetStat(tmpHetLeft);
        tmpPi2->setRightHetStat(tmpHetRight);
      }

      else if(suffix == "C")
      {
        auto tmpHetLeft = std::dynamic_pointer_cast<HetMoment>(getHetMoment(p1, p2, "B"));
        auto tmpHetRight = std::dynamic_pointer_cast<HetMoment>(getHetMoment(p3, p4, "A"));
        tmpPi2->setLeftHetStat(tmpHetLeft);
        tmpPi2->setRightHetStat(tmpHetRight);
      }

      else if(suffix == "D")
      {
        auto tmpHetLeft = std::dynamic_pointer_cast<HetMoment>(getHetMoment(p1, p2, "B"));
        auto tmpHetRight = std::dynamic_pointer_cast<HetMoment>(getHetMoment(p3, p4, "B"));
        tmpPi2->setLeftHetStat(tmpHetLeft);
        tmpPi2->setRightHetStat(tmpHetRight);
      }
    }
  }
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
