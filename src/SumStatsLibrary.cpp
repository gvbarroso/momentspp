/*
 * Authors: Gustavo V. Barroso
 * Created: 05/08/2022
 * Last modified: 13/02/2023
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

void SumStatsLibrary::aliasMoments(const std::vector<size_t>& selectedPopIds) // we assume selection acts on the left locus
{
  for(size_t i = 0; i < numDDStats_; ++i)
  {
    assert(moments_[i]->getPrefix() == "DD");

    size_t pop1 = moments_[i]->getPopIndices()[0];
    size_t pop2 = moments_[i]->getPopIndices()[1];

    if(pop1 != pop2) // cross-pop DD stats are aliased independently of the selection model
      moments_[i]->insertAlias(getDdMoment(pop2, pop1));
  }

  for(size_t i = numDDStats_; i < (numDDStats_ + numDzStats_); ++i)
  {
    assert(moments_[i]->getPrefix() == "Dz");

    size_t pop1 = moments_[i]->getPopIndices()[0];
    size_t pop2 = moments_[i]->getPopIndices()[1];
    size_t pop3 = moments_[i]->getPopIndices()[2];

    if(std::find(std::begin(selectedPopIds), std::end(selectedPopIds), pop2) == std::end(selectedPopIds)) // left locus is NOT under selection in pop2
      moments_[i]->insertAlias(getDzMoment(pop1, pop3, pop2));
  }

  for(size_t i = (numDDStats_ + numDzStats_); i < (numDDStats_ + numDzStats_ + numHetStats_); ++i)
  {
    assert(moments_[i]->getPrefix() == "H");

    auto tmp = std::dynamic_pointer_cast<HetMoment>(moments_[i]);

    size_t pop1 = tmp->getPopIndices()[0]; // left locus, potentially under selection
    size_t pop2 = tmp->getPopIndices()[1]; // right locus, always neutral by construction of summary statistics stored in Data class
    std::string suffix = tmp->getSuffix();

    if(pop1 == pop2) // within-population H's are always aliased
    {
      if(suffix == "A")
        moments_[i]->insertAlias(getHetMoment(pop1, pop2, "B"));

      else if(suffix == "B")
        moments_[i]->insertAlias(getHetMoment(pop1, pop2, "A"));
    }

    else // cross-population H's
    {
      if(tmp->isConstrained())
      {
        if(std::find(std::begin(selectedPopIds), std::end(selectedPopIds), pop1) != std::end(selectedPopIds)) // locus is under selection in pop1
        {
          if(suffix == "A") // pop1 carries derived allele in this HetMoment
            moments_[i]->insertAlias(getHetMoment(pop2, pop1, "B"));
        }

        else if(std::find(std::begin(selectedPopIds), std::end(selectedPopIds), pop2) != std::end(selectedPopIds)) // locus is under selection in pop2
        {
          if(suffix == "B") // pop2 carries derived allele in this HetMoment
            moments_[i]->insertAlias(getHetMoment(pop2, pop1, "A"));
        }
      }

      else // the locus concerning this HetMoment is not putatively under selection, alias p(1-p) and (1-p)p
      {
        if(suffix == "A")
        {
          moments_[i]->insertAlias(getHetMoment(pop2, pop1, "A"));
          moments_[i]->insertAlias(getHetMoment(pop2, pop1, "B"));
          moments_[i]->insertAlias(getHetMoment(pop1, pop2, "B"));
        }

        else if(suffix == "B")
        {
          moments_[i]->insertAlias(getHetMoment(pop2, pop1, "B"));
          moments_[i]->insertAlias(getHetMoment(pop2, pop1, "A"));
          moments_[i]->insertAlias(getHetMoment(pop1, pop2, "A"));
        }
      }
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

        bool testLeft = 1; // are the left H stats the same between pi2(i) and pi2(j)?

        if(left1 != left2)
        {
          if(left1->getAliases().size() == left2->getAliases().size())
          {
            for(size_t k = 0; k < left1->getAliases().size(); ++k)
            {
              if(!left1->hasAlias(left2))
                testLeft = 0;

              else
              {
                for(size_t l = 0; l < left2->getAliases().size(); ++l)
                {
                  if(!left1->hasAlias(left2->getAliases()[k]) && left1 != left2->getAliases()[k])
                    testLeft = 0;
                }
              }
            }
          }

          else
            testLeft = 0;
        }

        bool testRight = 1; // are the right H stats the same between pi2(i) and pi2(j)?

        if(right1 != right2)
        {
          if(right1->getAliases().size() == right2->getAliases().size())
          {
            for(size_t k = 0; k < right1->getAliases().size(); ++k)
            {
              if(!right1->hasAlias(right2))
                testLeft = 0;

              else
              {
                for(size_t l = 0; l < right2->getAliases().size(); ++l)
                {
                  if(!right1->hasAlias(right2->getAliases()[k]) && right1 != right2->getAliases()[k])
                    testLeft = 0;
                }
              }
            }
          }

          else
            testRight = 0;
        }

        bool permutable = 0; // permutations share selective status
        if((left1->isConstrained() == right2->isConstrained()) && (left2->isConstrained() == right1->isConstrained()))
        {
          if(left1->hasSamePopIds(right2) && left2->hasSamePopIds(right1))
            permutable = 1;
        }

        if((testLeft && testRight) || permutable)
          moments_[i]->insertAlias(std::dynamic_pointer_cast<Pi2Moment>(moments_[j]));
      }
    }
  }
}

std::vector<std::shared_ptr<Moment>> SumStatsLibrary::fetchCompressedBasis()
{
  std::vector<std::shared_ptr<Moment>> ret(0);
  ret.reserve(moments_.size());

  for(size_t i = 0; i < moments_.size(); ++i)
  {
    bool hasUniqueExpectation = 1;

    for(size_t j = 0; j < ret.size(); ++j)
    {
      if(ret[j]->hasAlias(moments_[i]))
        hasUniqueExpectation = 0;
    }

    if(hasUniqueExpectation)
      ret.emplace_back(moments_[i]);
  }

  return ret;
}

void SumStatsLibrary::initMoments_(const std::map<size_t, std::shared_ptr<Population>>& popMap)
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

      // NOTE insert H moments with isPutativelySelected_ == true based on popMap
      // maybe give them suffixes C and D?
      // also, do we want to include the "special" H's C and D if they concern populations not in popMap of *this epoch?

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

  /* NOTE proposal:
   *
   * assemple std::vector<size_t> selectedPopIds from popMap
   * call aliasMoments(selectedPopIds); [also make aliasMoments private?]
   *
   *
   */
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
