/*
 * Authors: Gustavo V. Barroso
 * Created: 05/08/2022
 * Last modified: 29/08/2023
 *
 */

#include <ios>

#include "SumStatsLibrary.hpp"

std::shared_ptr<Moment> SumStatsLibrary::getMoment(const std::string& name) const
{
  auto ptr = std::find_if(std::begin(moments_), std::end(moments_), [=](std::shared_ptr<Moment> tmp) { return tmp->getName() == name; });

  if(ptr != std::end(moments_))
    return *ptr;

  else
    throw bpp::Exception("SumStatsLibrary::could not find stat " + name);

  /*std::shared_ptr<Moment> ptr = nullptr;
  for(auto itMom = std::begin(moments_); itMom != std::end(moments_); ++itMom)
  {
    if((*itMom)->getName() == name)
      ptr = *itMom;
  }

  assert(ptr != nullptr);
  return ptr;*/
}

std::shared_ptr<Moment> SumStatsLibrary::getMoment(const std::string& prefix, const std::vector<size_t>& popIds, const std::vector<size_t>& factorIds) const
{
  return getMoment(assembleName_(prefix, popIds, factorIds));
}

std::shared_ptr<Moment> SumStatsLibrary::getMoment(size_t pos) const
{
  return basis_[pos];
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

void SumStatsLibrary::dropFactorIds(std::vector<size_t>& factorIds, size_t focalPopId, int removeCount) const
{
  assert(std::count(std::begin(factorIds), std::end(factorIds), focalPopId) >= removeCount);

  int counter = 0;
  for(auto it = std::begin(factorIds); it != std::end(factorIds);)
  {
    if(*it == focalPopId && counter < removeCount)
    {
      it = factorIds.erase(it);
      ++counter;
    }

    else
      ++it;
  }
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

void SumStatsLibrary::initMoments_(bool compress)
{
  moments_.reserve(getNumStats());

  size_t comb = 0;
  for(size_t i = 1; i < factorOrder_; ++i)
    comb += i;

  factorComb_ = comb;

  for(auto itI = std::begin(popIndices_); itI != std::end(popIndices_); ++itI)
  {
    for(auto itJ = std::begin(popIndices_); itJ != std::end(popIndices_); ++itJ)
    {
      std::string name = "DD_" + asString(*itI) + "_" + asString(*itJ);
      moments_.emplace_back(std::make_shared<DdMoment>(name, 0.)); // from the canonical Ragsdale-Gravel basis

      for(size_t i = 1; i < (factorOrder_ + 1); ++i) // for each factor id count
      {
        std::vector<size_t> factorIds(i);

        for(size_t j = 0; j < popIndices_.size(); ++j) // for each pop
        {
          for(size_t k = 0; k < i; ++k) // init factor ids to focal pop id in all positions
            factorIds[k] = popIndices_[j];

          for(size_t k = 0; k < i; ++k) // for each position in factor ids vector
          {
            for(size_t l = 0; l < popIndices_.size(); ++l)
            {
              factorIds[k] = popIndices_[l]; // switches pop id at position

              std::string nome = name + "_l";
              for(size_t m = 0; m < i; ++m)
                nome = nome + "_" + asString(factorIds[m]);

              moments_.emplace_back(std::make_shared<DdMoment>(nome, 0.));
            }
          }
        }
      }

      name = "Dr_" + asString(*itI) + "_" + asString(*itJ); // D_i_(1-2q)_j, where q is the freq of derived (neutral) allele in the right locus
      moments_.emplace_back(std::make_shared<DrMoment>(name, 0.));

      for(size_t i = 1; i < (factorOrder_ + 2); ++i) // NOTE Dr stats include one factor of (1-2p) more than other stats
      {
        std::vector<size_t> factorIds(i);

        for(size_t j = 0; j < popIndices_.size(); ++j)
        {
          for(size_t k = 0; k < i; ++k)
            factorIds[k] = popIndices_[j];

          for(size_t k = 0; k < i; ++k)
          {
            for(size_t l = 0; l < popIndices_.size(); ++l)
            {
              factorIds[k] = popIndices_[l];

              std::string nome = name + "_l";
              for(size_t m = 0; m < i; ++m)
                nome = nome + "_" + asString(factorIds[m]);

              moments_.emplace_back(std::make_shared<DrMoment>(nome, 0.));
            }
          }
        }
      }

      name = "Hl_" + asString(*itI) + "_" + asString(*itJ);
      moments_.emplace_back(std::make_shared<HetMoment>(name, 0., true)); // Hl_01 = p_0(1-p_1); Hl_10 = p_1(1-p_0)

      for(size_t i = 1; i < (factorOrder_ + 1); ++i)
      {
        std::vector<size_t> factorIds(i);

        for(size_t j = 0; j < popIndices_.size(); ++j)
        {
          for(size_t k = 0; k < i; ++k)
            factorIds[k] = popIndices_[j];

          for(size_t k = 0; k < i; ++k)
          {
            for(size_t l = 0; l < popIndices_.size(); ++l)
            {
              factorIds[k] = popIndices_[l];

              std::string nome = name + "_l";
              for(size_t m = 0; m < i; ++m)
                nome = nome + "_" + asString(factorIds[m]);

              moments_.emplace_back(std::make_shared<HetMoment>(nome, 0., true));
            }
          }
        }
      }

      // Hr stats don't require training factors of (1-2p_x)
      name = "Hr_" + asString(*itI) + "_" + asString(*itJ);
      moments_.emplace_back(std::make_shared<HetMoment>(name, 0., false)); // Hr_01 = p_0(1-p_1); Hr_10 = p_1(1-p_0)

      for(auto itK = std::begin(popIndices_); itK != std::end(popIndices_); ++itK)
      {
        for(auto itL = std::begin(popIndices_); itL != std::end(popIndices_); ++itL)
        {
          name = "pi2_" + asString(*itI) + "_" + asString(*itJ) + "_" + asString(*itK) + "_" + asString(*itL);
          moments_.emplace_back(std::make_shared<Pi2Moment>(name, 0., nullptr, nullptr));

          for(size_t i = 1; i < (factorOrder_ + 1); ++i)
          {
            std::vector<size_t> factorIds(i);

            for(size_t j = 0; j < popIndices_.size(); ++j)
            {
              for(size_t k = 0; k < i; ++k)
                factorIds[k] = popIndices_[j];

              for(size_t k = 0; k < i; ++k)
              {
                for(size_t l = 0; l < popIndices_.size(); ++l)
                {
                  factorIds[k] = popIndices_[l];

                  std::string nome = name + "_l";
                  for(size_t m = 0; m < i; ++m)
                    nome = nome + "_" + asString(factorIds[m]);

                  moments_.emplace_back(std::make_shared<Pi2Moment>(nome, 0., nullptr, nullptr));
                }
              }
            }
          }
        }
      }
    }
  }

  // includes "Dummy" Moment to convert into a homogeneous system (see Mutation::setUpMatrices_())
  moments_.emplace_back(std::make_shared<Moment>("I", 1.));
  countMoments_();
  printMoments(std::cout);

  for(size_t i = 0; i < moments_.size(); ++i)
    moments_[i]->setPosition(i);

  linkPi2HetStats_();
  basis_ = moments_; // default

  if(compress)
  {
    aliasMoments_();
    compressBasis_();
  }
}

std::string SumStatsLibrary::assembleName_(const std::string& prefix, const std::vector<size_t>& popIds, const std::vector<size_t>& factorIds) const
{
  /*std::cout << prefix << ";";
  for(auto& i : popIds)
    std::cout << i << ",";

  std::cout << ";";
  for(auto& i : factorIds)
    std::cout << i << ",";*/

  std::string name = prefix;
  for(size_t i = 0; i < popIds.size(); ++i)
    name = name + "_" + asString(popIds[i]);

  if(factorIds.size() > 0)
  {
    for(size_t i = 0; i < popIndices_.size(); ++i)
    {
      size_t count = std::count(std::begin(factorIds), std::end(factorIds), popIndices_[i]);

      if(count > 0)
        name = name + "_(1-2p" + asString(popIndices_[i]) + ")^" + asString(count);
    }
  }

  std::cout << name << "\n";
  return name;
}

void SumStatsLibrary::countMoments_()
{
  assert(moments_.size() != 0);

  // deletes duplicates introduced by clumsy creation of moments
  // and determines the ascending lexicographical order of stats in the rows/cols of matrices inside AbstractOperators
  std::sort(std::begin(moments_), std::end(moments_), compareMoments_);
  moments_.erase(std::unique(std::begin(moments_), std::end(moments_),
                             [=](std::shared_ptr<Moment> a, std::shared_ptr<Moment> b)
                             { return (a->getName() == b->getName() ? true : false); }), std::end(moments_));

  numDDStats_ = std::count_if(std::begin(moments_), std::end(moments_), [=] (std::shared_ptr<Moment> m) { return m->getPrefix() == "DD"; });
  numDrStats_ = std::count_if(std::begin(moments_), std::end(moments_), [=] (std::shared_ptr<Moment> m) { return m->getPrefix() == "Dr"; });
  numPi2Stats_ = std::count_if(std::begin(moments_), std::end(moments_), [=] (std::shared_ptr<Moment> m) { return m->getPrefix() == "pi2"; });
  numHetLeftStats_ = std::count_if(std::begin(moments_), std::end(moments_), [=] (std::shared_ptr<Moment> m) { return m->getPrefix() == "Hl"; });
  numHetRightStats_ = std::count_if(std::begin(moments_), std::end(moments_), [=] (std::shared_ptr<Moment> m) { return m->getPrefix() == "Hr"; });

  ///std::cout << numDDStats_ << "," << numDrStats_  << "," << numPi2Stats_  << "," << numHetLeftStats_  << "," << numHetRightStats_  << "\n";
}

void SumStatsLibrary::linkPi2HetStats_()
{
  for(size_t i = (1 + numDDStats_ + numDrStats_ + numHetLeftStats_ + numHetRightStats_); i < getNumStats(); ++i)
  {
    auto tmp = std::dynamic_pointer_cast<Pi2Moment>(moments_[i]);
    assert(tmp != nullptr);

    std::vector<size_t> dummy(0); // for searching Hr stat
    std::vector<size_t> popsLeft(0);
    std::vector<size_t> popsRight(0);

    popsLeft.reserve(2);
    popsRight.reserve(2);

    popsLeft.emplace_back(tmp->getPopIndices()[0]);
    popsLeft.emplace_back(tmp->getPopIndices()[1]);
    popsRight.emplace_back(tmp->getPopIndices()[2]);
    popsRight.emplace_back(tmp->getPopIndices()[3]);

    auto tmpHetLeft = std::dynamic_pointer_cast<HetMoment>(getMoment("Hl", popsLeft, tmp->getFactorIndices()));
    auto tmpHetRight = std::dynamic_pointer_cast<HetMoment>(getMoment("Hr", popsRight, dummy));

    assert(tmpHetLeft != nullptr && tmpHetRight != nullptr);

    tmp->setLeftHetStat(tmpHetLeft);
    tmp->setRightHetStat(tmpHetRight);
  }
}

void SumStatsLibrary::aliasMoments_() // selection acts on the left locus by design
{
  assert((numDDStats_ > 0) && (numDrStats_ > 0) & (numHetLeftStats_ > 0) && (numHetRightStats_ > 0) && (numPi2Stats_ > 0));

  for(size_t i = 0; i < (numDDStats_ + numDrStats_ + numHetLeftStats_ + numHetRightStats_); ++i)
  {
    std::vector<size_t> pops(0);
    pops.reserve(2);

    // finding candidate alias moments by inverting order of pop ids, NOTE: does this work for Dr moments as well?
    pops.emplace_back(moments_[i]->getPopIndices()[1]);
    pops.emplace_back(moments_[i]->getPopIndices()[0]);

    auto candidate = getMoment(assembleName_(moments_[i]->getPrefix(), pops, moments_[i]->getFactorIndices()));

    if(pops[0] != pops[1])
      moments_[i]->insertAlias(candidate);
  }

  // Pi2 stats are aliased if both their left and right H's are aliased OR if there's a left-right permutation
  for(size_t i = (numDDStats_ + numDrStats_ + numHetLeftStats_ + numHetRightStats_ + 1); i < getNumStats(); ++i)
  {
    auto left1 = std::dynamic_pointer_cast<Pi2Moment>(moments_[i])->getLeftHetStat();
    auto right1 = std::dynamic_pointer_cast<Pi2Moment>(moments_[i])->getRightHetStat();

    assert(left1 != nullptr && right1 != nullptr);

    for(size_t j = (1 + numDDStats_ + numDrStats_ + numHetLeftStats_ + numHetRightStats_); j < getNumStats(); ++j)
    {
      if(i != j)
      {
        assert(moments_[j]->getPrefix() == "pi2");

        auto left2 = std::dynamic_pointer_cast<Pi2Moment>(moments_[j])->getLeftHetStat();
        auto right2 = std::dynamic_pointer_cast<Pi2Moment>(moments_[j])->getRightHetStat();

        assert(left2 != nullptr && right2 != nullptr);

        bool leftEq = left1 == left2 || left1->hasAlias(left2);
        bool rightEq = right1 == right2 || right1->hasAlias(right2);

        bool leftRightPermut = 0; // permutations share selective status
        // NOTE CANNOT be true ATM because all HetLeft stats are being treated as necessarily constrained
        if((left1->isConstrained() == right2->isConstrained()) && (left2->isConstrained() == right1->isConstrained()))
        {
          if(left1->hasAlias(right2) && left2->hasAlias(right1))
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
  // reset, default was basis_ = moments_
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

