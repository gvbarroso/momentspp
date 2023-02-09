/*
 * Authors: Gustavo V. Barroso
 * Created: 05/08/2022
 * Last modified: 09/02/2023
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

void SumStatsLibrary::aliasMoments(const std::vector<size_t>& selectedPopIds)
{
  // NOTE selection acts on the left locus and may be epoch-specific
  // DD stats are aliased independently of the selection model
  for(size_t i = 0; i < numDDStats_; ++i)
  {
    assert(moments_[i]->getPrefix() == "DD");

    size_t pop1 = moments_[i]->getPopIndices()[0];
    size_t pop2 = moments_[i]->getPopIndices()[1];

    if(pop1 != pop2) // cross-pop D covar
      moments_[i]->insertAlias(getDdMoment(pop1, pop2)); // searching for alias, flip order
  }

  // Dz stats [D(1-2p)(1-2q)] are never aliased
  // H stats are aliased if
  for(size_t i = (numDDStats_ + numDzStats_); i < (numDDStats_ + numDzStats_ + numHetStats_); ++i)
  {
    assert(moments_[i]->getPrefix() == "H");

    auto tmp = std::dynamic_pointer_cast<HetMoment>(moments_[i]);

    // pop ids concerning allele frequencies:
    size_t pop1 = tmp->getPopIndices()[0]; // left locus, potentially under selection
    size_t pop2 = tmp->getPopIndices()[1]; // right locus, always neutral (by construction of summary statistics stored in Data class)
    std::string suffix = tmp->getSuffix();

    if(pop1 == pop2) // within-population H's always alias
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

        else // none of the populations represented has selection against the derived allele, alias p(1-p) and (1-p)p
        {
          if(suffix == "A")
            moments_[i]->insertAlias(getHetMoment(pop2, pop1, "B"));

          else if(suffix == "B")
            moments_[i]->insertAlias(getHetMoment(pop2, pop1, "A"));
        }
      }

      else // the locus concerning this HetMoment is not putatively under selection, alias p(1-p) and (1-p)p
      {
        if(suffix == "A")
          moments_[i]->insertAlias(getHetMoment(pop2, pop1, "B"));

        else if(suffix == "B")
          moments_[i]->insertAlias(getHetMoment(pop2, pop1, "A"));
      }
    }
  }

  // Pi2 stats are aliased if both their left and right H's are aliased
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

        bool testLeft = 1; // inits boolean to true,

        if(left1->getAliases().size() == left2->getAliases().size())
        {
          for(size_t k = 0; k < left1->getAliases().size(); ++k)
          {
            if(left1 != left2)
            {
              if(left1->getAliases()[k] != left2) // H's are potential aliases of each other
                testLeft = 0;
            }
          }
        }

        else
          testLeft = 0;

        bool testRight = 1; // inits boolean to true

        if(right1->getAliases().size() == right2->getAliases().size())
        {
          for(size_t k = 0; k < right1->getAliases().size(); ++k)
          {
            if(right1 != right2)
            {
              if(right1->getAliases()[k] != right2) // H's are potential aliases of each other
                testRight = 0;
            }
          }
        }

        else
          testRight = 0;

        if(testLeft && testRight)
          moments_[i]->insertAlias(std::dynamic_pointer_cast<Pi2Moment>(moments_[j]));
      }
    }
  }
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
