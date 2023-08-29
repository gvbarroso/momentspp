/*
 * Authors: Gustavo V. Barroso
 * Created: 05/08/2022
 * Last modified: 29/08/2023
 */


#ifndef _SUMSTATSLIBRARY_H_
#define _SUMSTATSLIBRARY_H_

#include <set>
#include <string>
#include <vector>
#include <memory>
#include <algorithm>
#include <cstring>
#include <utility>
#include <ostream>
#include <cassert>

#include <eigen3/Eigen/Core>

#include <Bpp/Text/TextTools.h>

#include "Population.hpp"
#include "Moment.hpp"
#include "DdMoment.hpp"
#include "DrMoment.hpp"
#include "HetMoment.hpp"
#include "Pi2Moment.hpp"


// intent is to have one instance of SumStatsLibrary per Epoch because each Epoch potentially has unique population sets
class SumStatsLibrary
{

private:
  size_t numPops_;
  size_t numDDStats_;
  size_t numDrStats_;
  size_t numHetLeftStats_; // left locus, potentially under selection
  size_t numHetRightStats_; // right locus, neutral
  size_t numPi2Stats_;
  size_t factorOrder_; // maximum number of 1-2p_x factors attached to a Moment
  size_t factorComb_; // number of ways we can attach 1-2p_x factors to a given moment, used to find moment indices quickly (NOTE working for P == 2)

  std::vector<size_t> popIndices_; // among all Moments in the Epoch to which *this belongs, stored for bookkeeping
  std::vector<std::shared_ptr<Moment>> moments_; // sorted alphabetically based on prefix_ and numerically based on popIndices_
  std::vector<std::shared_ptr<Moment>> basis_; // reduced # of moments, based on symmetries

public:
  SumStatsLibrary():
  numPops_(0),
  numDDStats_(0),
  numDrStats_(0),
  numHetLeftStats_(0),
  numHetRightStats_(0),
  numPi2Stats_(0),
  factorOrder_(0),
  factorComb_(0),
  popIndices_(0),
  moments_(0),
  basis_(0)
  { }

  SumStatsLibrary(const std::vector<std::shared_ptr<Population>>& pops, size_t factorOrder, bool compress):
  numPops_(pops.size()),
  numDDStats_(0),
  numDrStats_(0),
  numHetLeftStats_(0),
  numHetRightStats_(0),
  numPi2Stats_(0),
  factorOrder_(factorOrder),
  factorComb_(0),
  popIndices_(0),
  moments_(0),
  basis_(0)
  {
    popIndices_.reserve(pops.size());
    std::set<size_t> uniqueIds;

    for(auto it = std::begin(pops); it != std::end(pops); ++it)
    {
      popIndices_.emplace_back((*it)->getId());
      uniqueIds.insert((*it)->getId());
    }

    if(popIndices_.size() != uniqueIds.size())
      throw bpp::Exception("SumStatsLibrary::non-unique population indices!");

    std::sort(std::begin(popIndices_), std::end(popIndices_));
    initMoments_(compress);
  }

public:
  std::string asString(int i)
  {
    return bpp::TextTools::toString(i);
  }

  std::string asString(int i) const
  {
    return bpp::TextTools::toString(i);
  }

  size_t getNumPops()
  {
    return numPops_;
  }

  const std::vector<size_t>& getPopIndices() const
  {
    return popIndices_;
  }

  const std::vector<std::shared_ptr<Moment>>& getMoments() const
  {
    return moments_;
  }

  std::vector<std::shared_ptr<Moment>>& getMoments()
  {
    return moments_;
  }

  std::vector<std::shared_ptr<Moment>>& getBasis()
  {
    return basis_;
  }

  const std::vector<std::shared_ptr<Moment>>& getBasis() const
  {
    return basis_;
  }

  size_t getNumDDStats() const
  {
    return numDDStats_;
  }

  size_t getNumDrStats() const
  {
    return numDrStats_;
  }

  size_t getNumHetLeftStats() const
  {
    return numHetLeftStats_;
  }

  size_t getNumHetRightStats() const
  {
    return numHetRightStats_;
  }

  size_t getNumPi2Stats() const
  {
    return numPi2Stats_;
  }

  size_t getNumStats() const
  {
    return 1 + numDDStats_ + numDrStats_ + numHetLeftStats_ + numHetRightStats_ + numPi2Stats_;
  }

  size_t getSizeOfBasis() const
  {
    return basis_.size();
  }

  int getFactorOrder() const
  {
    return static_cast<int>(factorOrder_);
  }

  std::shared_ptr<Moment> getMoment(const std::string& name) const;

  std::shared_ptr<Moment> getMoment(const std::string& prefix, const std::vector<size_t>& popIds, const std::vector<size_t>& factorIds) const;

  std::shared_ptr<Moment> getMoment(size_t pos) const;

  std::shared_ptr<Moment> getDummyMoment() const
  {
    return moments_[getDummyIndexUncompressed()];
  }

  std::shared_ptr<Moment> getDummyMomentCompressed() const
  {
    return basis_[findCompressedIndex(getDummyIndexUncompressed())];
  }

  size_t getDummyIndexUncompressed() const
  {
    return numDDStats_ + numDrStats_ + numHetLeftStats_ + numHetRightStats_;
  }

  size_t findCompressedIndex(std::shared_ptr<Moment> mom) const;

  size_t findCompressedIndex(size_t uncompressedIndex) const;

  void dropFactorIds(std::vector<size_t>& factorIds, size_t focalPopId, int removeCount) const;

  Eigen::VectorXd fetchYvec();

  void printMoments(std::ostream& stream);

  void printBasis(std::ostream& stream);

private:
  void initMoments_(bool compress);

  static bool compareMoments_(std::shared_ptr<Moment> a, std::shared_ptr<Moment> b)
  {
    bool lessThan = 0;

    if(a->getPrefix() != b->getPrefix())
      lessThan = a->getPrefix() < b->getPrefix();

    else
    {
      auto aPops = a->getPopIndices();
      auto bPops = b->getPopIndices();

      assert(aPops.size() == bPops.size());

      if(aPops != bPops)
      {
        for(size_t i = 0; i < aPops.size(); ++i)
        {
          if(aPops[i] != bPops[i])
          {
            lessThan = aPops[i] < bPops[i];
            break;
          }
        }
      }

      else
      {
        if(a->getFactorPower() != b->getFactorPower())
          lessThan = a->getFactorPower() < b->getFactorPower();

        else
        {
          // these are sorted by design inside each Moment object
          auto aFactors = a->getFactorIndices();
          auto bFactors = b->getFactorIndices();

          assert(aFactors.size() == bFactors.size());

          if(aFactors != bFactors)
          {
            for(size_t i = 0; i < aFactors.size(); ++i)
            {
              if(aFactors[i] != bFactors[i])
              {
                lessThan = aFactors[i] < bFactors[i];
                break;
              }
            }
          }
        }
      }
    }

    return lessThan;
  }

  // for searching / comparing
  std::string assembleName_(const std::string& prefix, const std::vector<size_t>& popIds, const std::vector<size_t>& factorIds) const;

  void countMoments_();

  // assigns two HetMoment pointers to each Pi2Moment (left and right loci)
  void linkPi2HetStats_();

  // exploits symmetry among statistics to reduce size of basis, given constraints imposed by selection
  void aliasMoments_();

  void compressBasis_();

};

#endif

