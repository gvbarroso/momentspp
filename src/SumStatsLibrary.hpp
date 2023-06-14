/*
 * Authors: Gustavo V. Barroso
 * Created: 05/08/2022
 * Last modified: 14/06/2023
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
  size_t factorOrder_; // maximum number of 1-2p factors attached to a Moment

  std::vector<size_t> popIndices_; // among all Moments, stored for bookkeeping
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
  popIndices_(0),
  moments_(0),
  basis_(0)
  { }

  SumStatsLibrary(const std::vector<std::shared_ptr<Population>>& pops, size_t factorOrder, bool compressMoments):
  numPops_(pops.size()),
  numDDStats_(0),
  numDrStats_(0),
  numHetLeftStats_(0),
  numHetRightStats_(0),
  numPi2Stats_(0),
  factorOrder_(factorOrder),
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
    initMoments_(pops, compressMoments);
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

  size_t getFactorOrder() const
  {
    return factorOrder_;
  }

  std::shared_ptr<Moment> getMoment(const std::string& name) const;

  std::shared_ptr<Moment> getMoment(size_t pos) const;

  std::shared_ptr<DdMoment> getDdMoment(size_t id1, size_t id2, size_t factorPower) const;

  std::shared_ptr<DrMoment> getDrMoment(size_t id1, size_t id2, size_t factorPower) const;

  std::shared_ptr<HetMoment> getHetLeftMoment(size_t id1, size_t id2, size_t factorPower) const;

  std::shared_ptr<HetMoment> getHetRightMoment(size_t id1, size_t id2) const;

  std::shared_ptr<Pi2Moment> getPi2Moment(size_t id1, size_t id2, size_t id3, size_t id4, size_t factorPower) const;

  std::shared_ptr<Pi2Moment> getPi2Moment(std::shared_ptr<HetMoment> left, std::shared_ptr<HetMoment> right) const;

  std::shared_ptr<Moment> getDummyMoment() const
  {
    return moments_[getDummyIndexUncompressed()];
  }

  std::shared_ptr<Moment> getDummyMomentCompressed() const
  {
    return basis_[findCompressedIndex(getDummyIndexUncompressed())];
  }

  size_t findPopIndexRank(size_t index) const;

  size_t findDdIndex(size_t id1, size_t id2, size_t factorPower) const;

  size_t findDrIndex(size_t id1, size_t id2, size_t factorPower) const;

  size_t findHetLeftIndex(size_t id1, size_t id2, size_t factorPower) const;

  size_t findHetRightIndex(size_t id1, size_t id2) const;

  size_t getDummyIndexUncompressed() const
  {
    return numDDStats_ + numDrStats_ + numHetLeftStats_ + numHetRightStats_;
  }

  size_t findPi2Index(size_t id1, size_t id2, size_t id3, size_t id4, size_t factorPower) const;

  size_t findCompressedIndex(std::shared_ptr<Moment> mom) const;

  size_t findCompressedIndex(size_t uncompressedIndex) const;

  Eigen::VectorXd fetchYvec();

  void printMoments(std::ostream& stream);

  void printBasis(std::ostream& stream);

private:
  void initMoments_(const std::vector<std::shared_ptr<Population>>& pops, bool compressMoments);

  static bool compareMoments_(std::shared_ptr<Moment> a, std::shared_ptr<Moment> b)
  {
    bool lessThan = 0;

    if(a->getPrefix() != b->getPrefix())
      lessThan = a->getPrefix() < b->getPrefix();

    else
    {
      auto x = a->getPopIndices();
      auto y = b->getPopIndices();

      assert(x.size() == y.size());

      if(x != y)
      {
        for(size_t i = 0; i < x.size(); ++i)
        {
          if(x[i] != y[i])
          {
            lessThan = x[i] < y[i];
            break;
          }
        }
      }

      else
        lessThan = a->getFactorPower() < b->getFactorPower();
    }

    return lessThan;
  }

  void countMoments_();

  // assign two HetMoment pointers to each Pi2Moment (left and right loci)
  void linkPi2HetStats_();

  // the next two method exploit symmetry among statistics to reduce dimension of stats_, given constraints imposed by Selection:
  void aliasMoments_(const std::vector<size_t>& selectedPopIds);

  void compressBasis_();

};

#endif

