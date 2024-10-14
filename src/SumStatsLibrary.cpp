/*
 * Authors: Gustavo V. Barroso
 * Created: 05/08/2022
 * Last modified: 03/06/2024
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
  size_t ret = getSizeOfBasis(); // init to out-of-bounds

  for(size_t j = 0; j < basis_.size(); ++j)
  {
    if(basis_[j]->getName() == mom->getName() || basis_[j]->hasAlias(mom))
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

Eigen::Matrix<long double, Eigen::Dynamic, 1> SumStatsLibrary::fetchYvec()
{
  Eigen::Matrix<long double, Eigen::Dynamic, 1> y(moments_.size());

  for(size_t i = 0; i < moments_.size(); ++i)
    y(i) = moments_[i]->getValue();

  return y;
}

void SumStatsLibrary::printMoments(std::ostream& stream)
{
  for(size_t i = 0; i < moments_.size(); ++i)
  {
    stream << moments_[i] << "\t";
    moments_[i]->printAttributes(stream);
  }
}

void SumStatsLibrary::printBasis(std::ostream& stream)
{
  for(size_t i = 0; i < basis_.size(); ++i)
  {
    stream << basis_[i] << "\t";
    basis_[i]->printAttributes(stream);
  }
}

void SumStatsLibrary::readStatsFromFile(const std::string& fileName)
{
  std::cout << "\nReading statistics from file " << fileName << "\n";
  assert(moments_.size() > 0 && basis_.size() > 0);

  boost::iostreams::filtering_istream boostStream;
  std::ifstream inFile(fileName, std::ios_base::in);

  if(!inFile.is_open())
    throw bpp::Exception("SumStatsLibrary::could not open stats. file: " + fileName);

  else
    boostStream.push(inFile);

  std::string line = "";
  std::vector<std::string> splitLine(0);

  while(std::getline(boostStream, line))
  {
    // assumes "mom = val" are space-separated
    boost::split(splitLine, line, [](char c) { return c == ' '; });

    std::string name = splitLine[0];
    long double val = std::stod(splitLine[2]);

    auto mom = getMoment(name);
    mom->setValue(val);
  }
}

void SumStatsLibrary::initMoments_(bool compress)
{
  moments_.reserve(getNumStats());

  for(auto itI = std::begin(popIndices_); itI != std::end(popIndices_); ++itI)
  {
    #ifdef NAKED_D
    std::string nameD = "D_" + asString(*itI); // naked singed D
    moments_.emplace_back(std::make_shared<Moment>(nameD, 0.));

    for(size_t i = 1; i < (factorOrder_ + 3); ++i) // NOTE D stats include two factors of (1-2p) more than other stats
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

            std::string nome = nameD + "_l";
            for(size_t m = 0; m < i; ++m)
              nome = nome + "_" + asString(factorIds[m]);

            moments_.emplace_back(std::make_shared<Moment>(nome, 0.));
          }
        }
      }
    }
    #endif

    for(auto itJ = std::begin(popIndices_); itJ != std::end(popIndices_); ++itJ)
    {
      std::string name = "DD_" + asString(*itI) + "_" + asString(*itJ);
      moments_.emplace_back(std::make_shared<DdMoment>(name, 0.));

      for(size_t i = 1; i < (factorOrder_ + 1); ++i)
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

      for(size_t i = 1; i < (factorOrder_ + 3); ++i) // NOTE Dr stats include two factors of (1-2p) more than other stats
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

  cleanBasis_();
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

  return name;
}

void SumStatsLibrary::cleanBasis_()
{
  assert(moments_.size() != 0);

  // determines the ascending lexicographical order of stats in the rows/cols of matrices inside AbstractOperators
  std::sort(std::begin(moments_), std::end(moments_), compareMoments_);

  // deletes duplicates introduced by clumsy creation of moments
  moments_.erase(std::unique(std::begin(moments_), std::end(moments_),
                             [=](std::shared_ptr<Moment> a, std::shared_ptr<Moment> b)
                             { return (a->getName() == b->getName() ? true : false); }), std::end(moments_));

  for(size_t i = 0; i < moments_.size(); ++i)
    moments_[i]->setPosition(i);
}

void SumStatsLibrary::linkPi2HetStats_()
{
  for(size_t i = 0; i < getNumStats(); ++i)
  {
    auto tmpPi2 = std::dynamic_pointer_cast<Pi2Moment>(moments_[i]);

    if(tmpPi2 != nullptr)
    {
      std::vector<size_t> popsLeft(0);
      std::vector<size_t> popsRight(0);

      popsLeft.reserve(2);
      popsRight.reserve(2);

      popsLeft.emplace_back(tmpPi2->getPopIndices()[0]);
      popsLeft.emplace_back(tmpPi2->getPopIndices()[1]);
      popsRight.emplace_back(tmpPi2->getPopIndices()[2]);
      popsRight.emplace_back(tmpPi2->getPopIndices()[3]);

      auto tmpHetLeft = std::dynamic_pointer_cast<HetMoment>(getMoment("Hl", popsLeft, tmpPi2->getFactorIndices()));
      auto tmpHetRight = std::dynamic_pointer_cast<HetMoment>(getMoment("Hr", popsRight, { }));

      assert(tmpHetLeft != nullptr && tmpHetRight != nullptr);

      tmpPi2->setLeftHetStat(tmpHetLeft);
      tmpPi2->setRightHetStat(tmpHetRight);
    }
  }
}

void SumStatsLibrary::aliasMoments_() // selection acts on the left locus by design
{
  assert(getNumStats() > 0);

  for(size_t i = 0; i < getNumStats(); ++i) // NOTE alias Hl_ij IFF s_i == s_j
  {
    if((moments_[i]->getPrefix() == "DD") || (moments_[i]->getPrefix() == "Hr"))
    {
      std::vector<size_t> pops(0);
      pops.reserve(2);

      // finding candidate alias moments by inverting order of pop ids
      pops.emplace_back(moments_[i]->getPopIndices()[1]);
      pops.emplace_back(moments_[i]->getPopIndices()[0]);

      auto candidate = getMoment(assembleName_(moments_[i]->getPrefix(), pops, moments_[i]->getFactorIndices()));

      if(pops[0] != pops[1])
        moments_[i]->insertAlias(candidate);
    }
  }

  // pi2 stats are aliased if their Hl are the same and and their Hr are aliased
  for(size_t i = 0; i < getNumStats(); ++i)
  {
    if(moments_[i]->getPrefix() == "pi2")
    {
      auto tmpPi2First = std::dynamic_pointer_cast<Pi2Moment>(moments_[i]);

      auto left1 = tmpPi2First->getLeftHetStat();
      auto right1 = tmpPi2First->getRightHetStat();

      assert(left1 != nullptr && right1 != nullptr);

      for(size_t j = i + 1; j < getNumStats(); ++j)
      {
        if(moments_[j]->getPrefix() == "pi2")
        {
          auto tmpPi2Second = std::dynamic_pointer_cast<Pi2Moment>(moments_[j]);

          auto left2 = tmpPi2Second->getLeftHetStat();
          auto right2 = tmpPi2Second->getRightHetStat();

          assert(left2 != nullptr && right2 != nullptr);

          bool leftEq = left1 == left2 || left1->hasAlias(left2);
          bool rightEq = right1 == right2 || right1->hasAlias(right2);

          if(leftEq && rightEq)
            tmpPi2First->insertAlias(tmpPi2Second);
        }
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
