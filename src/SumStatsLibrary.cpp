/*
 * Authors: Gustavo V. Barroso
 * Created: 05/08/2022
 * Last modified: 14/06/2023
 *
 */

#include <ios>

#include "SumStatsLibrary.hpp"

std::shared_ptr<Moment> SumStatsLibrary::getMoment(const std::string& name) const
{
  std::shared_ptr<Moment> ptr = nullptr;
  for(auto itMom = std::begin(moments_); itMom != std::end(moments_); ++itMom)
  {
    if((*itMom)->getName() == name)
      ptr = *itMom;
  }

  assert(ptr != nullptr);
  return ptr;
}

std::shared_ptr<Moment> SumStatsLibrary::getMoment(size_t pos) const
{
  return basis_[pos];
}

std::shared_ptr<DdMoment> SumStatsLibrary::getDdMoment(size_t id1, size_t id2, size_t factorPower) const
{
  size_t focalMomIndex = findDdIndex(id1, id2, factorPower);
  auto ret = std::dynamic_pointer_cast<DdMoment>(moments_[focalMomIndex]);

  if(ret != nullptr)
    return ret;

  else
    throw bpp::Exception("SumStatsLibrary::bad dynamic_pointer_cast attempt: DD" + asString(id1) + asString(id2));
}

std::shared_ptr<DrMoment> SumStatsLibrary::getDrMoment(size_t id1, size_t id2, size_t factorPower) const
{
  size_t focalMomIndex = findDrIndex(id1, id2, factorPower);
  auto ret = std::dynamic_pointer_cast<DrMoment>(moments_[focalMomIndex]);

  if(ret != nullptr)
    return ret;

  else
    throw bpp::Exception("SumStatsLibrary::bad dynamic_pointer_cast attempt: Dr" + asString(id1) + asString(id2));
}

std::shared_ptr<HetMoment> SumStatsLibrary::getHetLeftMoment(size_t id1, size_t id2, size_t factorPower) const
{
  size_t focalMomIndex = findHetLeftIndex(id1, id2, factorPower);
  auto ret = std::dynamic_pointer_cast<HetMoment>(moments_[focalMomIndex]);

  if(ret != nullptr)
    return ret;

  else
    throw bpp::Exception("SumStatsLibrary::bad dynamic_pointer_cast attempt: Hl" + asString(id1) + asString(id2));
}

std::shared_ptr<HetMoment> SumStatsLibrary::getHetRightMoment(size_t id1, size_t id2) const
{
  size_t focalMomIndex = findHetRightIndex(id1, id2);
  auto ret = std::dynamic_pointer_cast<HetMoment>(moments_[focalMomIndex]);

  if(ret != nullptr)
    return ret;

  else
    throw bpp::Exception("SumStatsLibrary::bad dynamic_pointer_cast attempt: Hr" + asString(id1) + asString(id2));
}

std::shared_ptr<Pi2Moment> SumStatsLibrary::getPi2Moment(size_t id1, size_t id2, size_t id3, size_t id4, size_t factorPower) const
{
  size_t focalMomIndex = findPi2Index(id1, id2, id3, id4, factorPower);
  auto ret = std::dynamic_pointer_cast<Pi2Moment>(moments_[focalMomIndex]);

  if(ret != nullptr)
    return ret;

  else
    throw bpp::Exception("SumStatsLibrary::bad dynamic_pointer_cast attempt: Pi2" + asString(id1) + asString(id2) + asString(id3) + asString(id3));
}

std::shared_ptr<Pi2Moment> SumStatsLibrary::getPi2Moment(std::shared_ptr<HetMoment> left, std::shared_ptr<HetMoment> right) const
{
  assert(left != nullptr && right != nullptr);
  std::shared_ptr<Pi2Moment> ret = nullptr;

  for(size_t i = (numDDStats_ + numDrStats_ + numHetLeftStats_ + numHetRightStats_ + 1); i < getNumStats(); ++i)
  {
    auto tmp = std::dynamic_pointer_cast<Pi2Moment>(moments_[i]);
    assert(tmp != nullptr);

    if(left == tmp->getLeftHetStat() && right == tmp->getRightHetStat()) // WARNING
      ret = tmp;
  }

  if(ret != nullptr)
    return ret;

  else
    throw bpp::Exception("SumStatsLibrary::could not find pi2 stat for " + left->getName() + " + " + right->getName());
}

size_t SumStatsLibrary::findPopIndexRank(size_t index) const // among all pop indices
{
  auto it = std::find(std::begin(popIndices_), std::end(popIndices_), index);
  assert(it != std::end(popIndices_));

  return std::distance(std::begin(popIndices_), it); // indexed from 0
}

size_t SumStatsLibrary::findDdIndex(size_t id1, size_t id2, size_t factorPower) const
{
  size_t r1 = findPopIndexRank(id1);
  size_t r2 = findPopIndexRank(id2);

  return r1 * numPops_ * (factorOrder_ + 1) + r2 * (factorOrder_ + 1) + factorPower;
}

size_t SumStatsLibrary::findDrIndex(size_t id1, size_t id2, size_t factorPower) const
{
  size_t r1 = findPopIndexRank(id1);
  size_t r2 = findPopIndexRank(id2);

  return numDDStats_ + r1 * numPops_ * (factorOrder_ + 1) + r2 * (factorOrder_ + 1) + factorPower;
}

size_t SumStatsLibrary::findHetLeftIndex(size_t id1, size_t id2, size_t factorPower) const
{
  size_t r1 = findPopIndexRank(id1);
  size_t r2 = findPopIndexRank(id2);

  return numDDStats_ + numDrStats_ + r1 * numPops_ * (factorOrder_ + 1) + r2 * (factorOrder_ + 1) + factorPower;
}

size_t SumStatsLibrary::findHetRightIndex(size_t id1, size_t id2) const
{
  size_t r1 = findPopIndexRank(id1);
  size_t r2 = findPopIndexRank(id2);

  return numDDStats_ + numDrStats_ + numHetLeftStats_ + r1 * numPops_ * (factorOrder_ + 1) + r2 * (factorOrder_ + 1);
}

size_t SumStatsLibrary::findPi2Index(size_t id1, size_t id2, size_t id3, size_t id4, size_t factorPower) const
{
  size_t r1 = findPopIndexRank(id1);
  size_t r2 = findPopIndexRank(id2);
  size_t r3 = findPopIndexRank(id3);
  size_t r4 = findPopIndexRank(id4);

  // 1 + because of dummy Moment "I_" after "H_**" to make system homogeneous(see initMoments_())
  return 1 + numDDStats_ + numDrStats_ + numHetLeftStats_ + numHetRightStats_ + r1 * numPops_ * numPops_ * numPops_ * (factorOrder_ + 1) + r2 * numPops_ * numPops_ * (factorOrder_ + 1) + r3 * numPops_ * (factorOrder_ + 1) + r4 * (factorOrder_ + 1) + factorPower;
}

size_t SumStatsLibrary::findCompressedIndex(std::shared_ptr<Moment> mom) const
{
  size_t ret = getNumStats(); // init to out-of-bounds

  for(size_t j = 0; j < basis_.size(); ++j)
  {
    if(basis_[j] == mom || basis_[j]->hasAlias(mom))
      ret = j;
  }

  if(ret < getNumStats())
    return ret;

  else
    throw bpp::Exception("SumStatsLibrary::could not find compressed index for " + mom->getName());
}

size_t SumStatsLibrary::findCompressedIndex(size_t uncompressedIndex) const
{
  size_t ret = getNumStats(); // init to out-of-bounds
  auto mom = moments_[uncompressedIndex];

  for(size_t j = 0; j < basis_.size(); ++j)
  {
    if(basis_[j] == mom || basis_[j]->hasAlias(mom))
      ret = j;
  }

  if(ret < getNumStats())
    return ret;

  else
    throw bpp::Exception("SumStatsLibrary::could not find compressed index for " + mom->getName());
}

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

void SumStatsLibrary::printBasis(std::ostream& stream)
{
  for(size_t i = 0; i < basis_.size(); ++i)
    basis_[i]->printAttributes(stream);
}

void SumStatsLibrary::initMoments_(const std::vector<std::shared_ptr<Population>>& pops, bool compress)
{
  moments_.reserve(getNumStats());

  // fetches ids of populations where derived allele for left locus is under selection
  std::vector<size_t> selectedPopIds(0);
  selectedPopIds.reserve(pops.size());
  for(auto it = std::begin(pops); it != std::end(pops); ++it)
  {
    if((*it)->hasSelection())
      selectedPopIds.emplace_back((*it)->getId());
  }

  for(auto itI = std::begin(popIndices_); itI != std::end(popIndices_); ++itI)
  {
    for(auto itJ = std::begin(popIndices_); itJ != std::end(popIndices_); ++itJ)
    {
      std::string name = "DD_" + asString(*itI) + "_" + asString(*itJ) + "_l";

      for(size_t i = 0; i < factorOrder_ + 1; ++i)
      {
        moments_.emplace_back(std::make_shared<DdMoment>(name, 0.));
        name = name + "_" + asString(selectedPopIds.front());
      }

      name = "Hl_" + asString(*itI) + "_" + asString(*itJ) + "_l";

      for(size_t i = 0; i < factorOrder_ + 1; ++i)
      {
        moments_.emplace_back(std::make_shared<HetMoment>(name, 0., true)); // Hl_01 = p_0(1-p_1); Hl_10 = p_1(1-p_0)
        name = name + "_" + asString(selectedPopIds.front());
      }

      name = "Hr_" + asString(*itI) + "_" + asString(*itJ);
      moments_.emplace_back(std::make_shared<HetMoment>(name, 0., false)); // Hr_01 = p_0(1-p_1); Hr_10 = p_1(1-p_0)

      name = "Dr_" + asString(*itI) + "_" + asString(*itJ) + "_l"; // D_i_(1-2q)_j, where q is the freq of derived (neutral) allele in the right locus

      for(size_t i = 0; i < factorOrder_ + 1; ++i)
      {
        moments_.emplace_back(std::make_shared<DrMoment>(name, 0.));
        name = name + "_" + asString(selectedPopIds.front());
      }

      for(auto itK = std::begin(popIndices_); itK != std::end(popIndices_); ++itK)
      {
        for(auto itL = std::begin(popIndices_); itL != std::end(popIndices_); ++itL)
        {
          name = "pi2_" + asString(*itI) + "_" + asString(*itJ) + "_" + asString(*itK) + "_" + asString(*itL) + "_l";

          for(size_t i = 0; i < factorOrder_ + 1; ++i)
          {
            moments_.emplace_back(std::make_shared<Pi2Moment>(name, 0., nullptr, nullptr));
            name = name + "_" + asString(selectedPopIds.front());
          }
        }
      }
    }
  }

  // includes "Dummy" Moment to convert into a homogeneous system (see Mutation::setUpMatrices_())
  moments_.emplace_back(std::make_shared<Moment>("I", 1.));

  // determines the ascending lexicographical order of stats in the rows/cols of matrices inside AbstractOperators
  std::sort(std::begin(moments_), std::end(moments_), compareMoments_);
  countMoments_();

  printMoments(std::cout);

  for(size_t i = 0; i < moments_.size(); ++i)
    moments_[i]->setPosition(i);

  linkPi2HetStats_();
  basis_ = moments_; // default

  /* NOTE activate when aliasMoments_ is fixed for P > 1
  if(compress)
  {
    aliasMoments_(selectedPopIds);
    compressBasis_();
  }
  */
}

void SumStatsLibrary::countMoments_()
{
  assert(moments_.size() != 0);

  numDDStats_ = std::count_if(std::begin(moments_), std::end(moments_), [] (std::shared_ptr<Moment> m) { return m->getPrefix() == "DD"; });
  numDrStats_ = std::count_if(std::begin(moments_), std::end(moments_), [] (std::shared_ptr<Moment> m) { return m->getPrefix() == "Dr"; });
  numPi2Stats_ = std::count_if(std::begin(moments_), std::end(moments_), [] (std::shared_ptr<Moment> m) { return m->getPrefix() == "pi2"; });
  numHetLeftStats_ = std::count_if(std::begin(moments_), std::end(moments_), [] (std::shared_ptr<Moment> m) { return m->getPrefix() == "Hl"; });
  numHetRightStats_ = std::count_if(std::begin(moments_), std::end(moments_), [] (std::shared_ptr<Moment> m) { return m->getPrefix() == "Hr"; });
}

void SumStatsLibrary::linkPi2HetStats_()
{
  for(size_t i = (1 + numDDStats_ + numDrStats_ + numHetLeftStats_ + numHetRightStats_); i < getNumStats(); ++i)
  {
    auto tmpPi2 = std::dynamic_pointer_cast<Pi2Moment>(moments_[i]);
    assert(tmpPi2 != nullptr);

    size_t p1 = tmpPi2->getPopIndices()[0];
    size_t p2 = tmpPi2->getPopIndices()[1];
    size_t p3 = tmpPi2->getPopIndices()[2];
    size_t p4 = tmpPi2->getPopIndices()[3];
    size_t x = tmpPi2->getFactorPower();

    auto tmpHetLeft = std::dynamic_pointer_cast<HetMoment>(getHetLeftMoment(p1, p2, x));
    auto tmpHetRight = std::dynamic_pointer_cast<HetMoment>(getHetRightMoment(p3, p4));

    assert(tmpHetLeft != nullptr && tmpHetRight != nullptr);

    tmpPi2->setLeftHetStat(tmpHetLeft);
    tmpPi2->setRightHetStat(tmpHetRight);
  }
}

/* NOTE needs to be fixed when P > 1
void SumStatsLibrary::aliasMoments_(const std::vector<size_t>& selectedPopIds) // selection acts on the left locus by design
{
  for(size_t i = 0; i < numDDStats_; ++i)
  {
    assert(moments_[i]->getPrefix() == "DD");

    size_t pop1 = moments_[i]->getPopIndices()[0];
    size_t pop2 = moments_[i]->getPopIndices()[1];

    bool statusP1 = std::find(std::begin(selectedPopIds), std::end(selectedPopIds), pop1) != std::end(selectedPopIds);
    bool statusP2 = std::find(std::begin(selectedPopIds), std::end(selectedPopIds), pop2) != std::end(selectedPopIds);

    if((pop1 != pop2) && (statusP1 == statusP2))
      moments_[i]->insertAlias(getDdMoment(pop2, pop1));
  }

  // TODO when P > 1, use factorIds for aliasing permutations with same selective constraint
  for(size_t i = numDDStats_; i < (numDDStats_ + numDrStats_); ++i)
  {
    assert(moments_[i]->getPrefix() == "Dr");

    size_t pop1 = moments_[i]->getPopIndices()[0];
    size_t pop2 = moments_[i]->getPopIndices()[1];

    auto factorIds = moments_[i]->getFactorIndices();

    for(size_t j = 0; j < factorIds.size(); ++j)
    {
      bool sel = std::find(std::begin(selectedPopIds), std::end(selectedPopIds), factorIds[j]) != std::end(selectedPopIds);

      if(sel)
      {
        moments_[i]->insertAlias(getDrMoment(pop1, p2));

    }
  }

  for(size_t i = (numDDStats_ + numDrStats_); i < (numDDStats_ + numDrStats_ + numHetLeftStats_ + numHetRightStats_); ++i)
  {
    assert(moments_[i]->getPrefix() == "H");

    auto tmp = std::dynamic_pointer_cast<HetMoment>(moments_[i]);

    size_t pop1 = tmp->getPopIndices()[0]; // derived allele, potentially under selection
    size_t pop2 = tmp->getPopIndices()[1]; // ancestral allele, neutral by construction

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
  for(size_t i = (numDDStats_ + numDrStats_ + numHetLeftStats_ + numHetRightStats_ + 1); i < getNumStats(); ++i)
  {
    assert(moments_[i]->getPrefix() == "pi2");

    auto left1 = std::dynamic_pointer_cast<Pi2Moment>(moments_[i])->getLeftHetStat();
    auto right1 = std::dynamic_pointer_cast<Pi2Moment>(moments_[i])->getRightHetStat();

    for(size_t j = (numDDStats_ + numDrStats_ + numHetLeftStats_ + numHetRightStats_ + 1); j < getNumStats(); ++j)
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
*/

void SumStatsLibrary::compressBasis_()
{
  basis_.clear();
  basis_.reserve(moments_.size());

  for(size_t i = 0; i < moments_.size(); ++i)
  {
    bool hasUniqueExpectation = 1;

    for(size_t j = 0; j < basis_.size(); ++j)
    {
      if(basis_[j]->hasAlias(moments_[i]))
        hasUniqueExpectation = 0;
    }

    if(hasUniqueExpectation)
      basis_.emplace_back(moments_[i]);
  }

  for(size_t i = 0; i < basis_.size(); ++i)
    basis_[i]->setPosition(i);
}

