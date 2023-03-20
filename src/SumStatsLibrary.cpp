/*
 * Authors: Gustavo V. Barroso
 * Created: 05/08/2022
 * Last modified: 20/03/2023
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
  for(size_t i = 0; i < compressedBasis_.size(); ++i)
    compressedBasis_[i]->printAttributes(stream);
}

void SumStatsLibrary::initMoments_(const std::map<size_t, std::shared_ptr<Population>>& popMap)
{
  moments_.reserve(getNumStats());

  // first pass to include DD, Dz, H statistics
  for(auto itI = std::begin(popIndices_); itI != std::end(popIndices_); ++itI)
  {
    for(auto itJ = std::begin(popIndices_); itJ != std::end(popIndices_); ++itJ)
    {
      moments_.emplace_back(std::make_shared<DdMoment>("DD_" + asString(*itI) + "_" + asString(*itJ), 0.));
      moments_.emplace_back(std::make_shared<HetMoment>("H_" + asString(*itI) + "_" + asString(*itJ), 0., false));
      // NOTE: H_01 = p_0(1-p_1); H_10 = p_1(1-p_0)
      // insert H moments with isPutativelySelected_ == true based on popMap

      for(auto itK = std::begin(popIndices_); itK != std::end(popIndices_); ++itK)
      {
        moments_.emplace_back(std::make_shared<DzMoment>("Dz_" + asString(*itI) + "_" + asString(*itJ) + "_" + asString(*itK), 0.));

        for(auto itL = std::begin(popIndices_); itL != std::end(popIndices_); ++itL)
          moments_.emplace_back(std::make_shared<Pi2Moment>("pi2_" + asString(*itI) + "_" + asString(*itJ) + "_" + asString(*itK) + "_" + asString(*itL), 0., nullptr, nullptr));
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

  // fetches ids of populations where derived allele for left locus is under selection
  std::vector<size_t> selectedPopIds(0);
  selectedPopIds.reserve(popMap.size());
  for(auto it = std::begin(popMap); it != std::end(popMap); ++it)
  {
    if(it->second->hasSelection())
      selectedPopIds.emplace_back(it->second->getId());
  }

  aliasMoments_(selectedPopIds);
  compressBasis_();
}

// for each Pi2Moment, sets the two pointers corresponding to HetMoments (left and right loci)
void SumStatsLibrary::linkPi2HetStats_()
{
  for(size_t i = (numDDStats_ + numDzStats_ + numHetStats_ + 1); i < getNumStats(); ++i)
  {
    auto tmpPi2 = std::dynamic_pointer_cast<Pi2Moment>(moments_[i]);
    assert(tmpPi2 != nullptr);

    size_t p1 = tmpPi2->getPopIndices()[0];
    size_t p2 = tmpPi2->getPopIndices()[1];
    size_t p3 = tmpPi2->getPopIndices()[2];
    size_t p4 = tmpPi2->getPopIndices()[3];

    auto tmpHetLeft = std::dynamic_pointer_cast<HetMoment>(getHetMoment(p1, p2));
    auto tmpHetRight = std::dynamic_pointer_cast<HetMoment>(getHetMoment(p3, p4));
    assert(tmpHetLeft != nullptr && tmpHetRight != nullptr);

    tmpPi2->setLeftHetStat(tmpHetLeft);
    tmpPi2->setRightHetStat(tmpHetRight);
  }
}

void SumStatsLibrary::aliasMoments_(const std::vector<size_t>& selectedPopIds) // we assume selection acts on the left locus
{
  for(size_t i = 0; i < numDDStats_; ++i)
  {
    assert(moments_[i]->getPrefix() == "DD");

    size_t pop1 = moments_[i]->getPopIndices()[0];
    size_t pop2 = moments_[i]->getPopIndices()[1];

    if(pop1 != pop2) // NOTE cross-pop DD stats are aliased independently of the selection model? check D^2 under selection, left and right
      moments_[i]->insertAlias(getDdMoment(pop2, pop1));
  }

  for(size_t i = numDDStats_; i < (numDDStats_ + numDzStats_); ++i)
  {
    assert(moments_[i]->getPrefix() == "Dz");

    size_t pop1 = moments_[i]->getPopIndices()[0];
    size_t pop2 = moments_[i]->getPopIndices()[1];
    size_t pop3 = moments_[i]->getPopIndices()[2];

    if(pop2 != pop3)
    {
      // if left locus is NOT under selection in pop2 (right locus is neutral by construction)
      if(std::find(std::begin(selectedPopIds), std::end(selectedPopIds), pop2) == std::end(selectedPopIds))
        moments_[i]->insertAlias(getDzMoment(pop1, pop3, pop2));
    }
  }

  for(size_t i = (numDDStats_ + numDzStats_); i < (numDDStats_ + numDzStats_ + numHetStats_); ++i)
  {
    assert(moments_[i]->getPrefix() == "H");

    auto tmp = std::dynamic_pointer_cast<HetMoment>(moments_[i]);

    size_t pop1 = tmp->getPopIndices()[0]; // left locus, potentially under selection
    size_t pop2 = tmp->getPopIndices()[1]; // right locus, always neutral by construction of summary statistics stored in Data class

    if(pop1 != pop2)
    {
      if(tmp->isConstrained()) // HetMoment concerns left locus, under selection in at least 1 pop
      {
        bool pop1Sel = std::find(std::begin(selectedPopIds), std::end(selectedPopIds), pop1) != std::end(selectedPopIds);
        bool pop2Sel = std::find(std::begin(selectedPopIds), std::end(selectedPopIds), pop1) != std::end(selectedPopIds);

        if(pop1Sel == pop2Sel) // if pop1 and pop2 share status of constraint
          moments_[i]->insertAlias(getHetMoment(pop2, pop1));
      }

      else // the locus concerning this HetMoment is not putatively under selection
        moments_[i]->insertAlias(getHetMoment(pop2, pop1));
    }
  }

  // Pi2 stats are aliased if both their left and right H's are aliased OR if focal pops don't experience selection and there's a left-right permutation
  for(size_t i = (numDDStats_ + numDzStats_ + numHetStats_ + 1); i < getNumStats(); ++i)
  {
    assert(moments_[i]->getPrefix() == "pi2");

    auto left1 = std::dynamic_pointer_cast<Pi2Moment>(moments_[i])->getLeftHetStat();
    auto right1 = std::dynamic_pointer_cast<Pi2Moment>(moments_[i])->getRightHetStat();

    for(size_t j = (numDDStats_ + numDzStats_ + numHetStats_ + 1); j < getNumStats(); ++j)
    {
      if(i != j)
      {
        assert(moments_[j]->getPrefix() == "pi2");

        auto left2 = std::dynamic_pointer_cast<Pi2Moment>(moments_[j])->getLeftHetStat();
        auto right2 = std::dynamic_pointer_cast<Pi2Moment>(moments_[j])->getRightHetStat();

        bool leftEq = left1 == left2 || left1->hasSamePopIds(left2);
        bool rightEq = right1 == right2 || right1->hasSamePopIds(right2);

        bool leftRightPermut = 0; // permutations share selective status
        if((left1->isConstrained() == right2->isConstrained()) && (left2->isConstrained() == right1->isConstrained()))
        {
          if(left1->hasSamePopIds(right2) && left2->hasSamePopIds(right1))
            leftRightPermut = 1;
        }

        if((leftEq && rightEq) || leftRightPermut)
          moments_[i]->insertAlias(std::dynamic_pointer_cast<Pi2Moment>(moments_[j]));
      }
    }
  }
}

void SumStatsLibrary::compressBasis_()
{
  compressedBasis_.reserve(moments_.size());

  for(size_t i = 0; i < moments_.size(); ++i)
  {
    bool hasUniqueExpectation = 1;

    for(size_t j = 0; j < compressedBasis_.size(); ++j)
    {
      if(compressedBasis_[j]->hasAlias(moments_[i]))
        hasUniqueExpectation = 0;
    }

    if(hasUniqueExpectation)
      compressedBasis_.emplace_back(moments_[i]);
  }

  for(size_t i = 0; i < compressedBasis_.size(); ++i)
    compressedBasis_[i]->setPosition(i);
}

