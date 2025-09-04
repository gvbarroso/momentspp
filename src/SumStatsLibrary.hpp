/*
 * Authors: Gustavo V. Barroso
 * Created: 05/08/2022
 * Last modified: 03/09/2025
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

#include <boost/algorithm/string.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/zlib.hpp>

#include <eigen3/Eigen/Core>
#include <eigen3/unsupported/Eigen/MPRealSupport> // for arbitrary-precision arithmetic

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
  size_t factorOrder_; // maximum number of 1-2p_x factors attached to a Moment

  std::vector<size_t> popIndices_; // among all Moments in the Epoch to which *this belongs, stored for bookkeeping
  std::vector<std::shared_ptr<Moment>> moments_; // sorted alphabetically based on prefix_ and numerically based on popIndices_
  std::vector<std::shared_ptr<Moment>> basis_; // reduced # of moments, compressed based on population genetic symmetries

public:
  SumStatsLibrary():
  numPops_(0),
  factorOrder_(0),
  popIndices_(0),
  moments_(0),
  basis_(0)
  { }

  SumStatsLibrary(const std::vector<std::shared_ptr<Population>>& pops, size_t factorOrder, bool compress):
  numPops_(pops.size()),
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

  size_t getSizeOfBasis() const
  {
    return basis_.size();
  }

  size_t getNumStats() const
  {
    return moments_.size();
  }

  size_t getFactorOrder() const
  {
    return factorOrder_;
  }

  std::shared_ptr<Moment> getMoment(const std::string& name) const;

  std::shared_ptr<Moment> getMoment(const std::string& prefix, const std::vector<size_t>& popIds, const std::vector<size_t>& factorIds) const;

  std::shared_ptr<Moment> getMoment(size_t pos) const;

  size_t findCompressedIndex(std::shared_ptr<Moment> mom) const;

  size_t findCompressedIndex(size_t uncompressedIndex) const;

  void dropFactorIds(std::vector<size_t>& factorIds, size_t focalPopId, int removeCount) const;

  Eigen::Matrix<mpfr::mpreal, Eigen::Dynamic, 1> fetchYvec();

  void printMoments(std::ostream& stream);

  void printBasis(std::ostream& stream);

  void readStatsFromFile(const std::string& fileName);

  // this assumes a maximum population count of 2 (speed constraint imposed by selection) and is used for convenience in the Drift operator
  size_t fetchOtherId(size_t id)
  {
    assert(popIndices_.size() == 2);

    if(popIndices_[0] == id)
      return popIndices_[1];

    else
      return popIndices_[0];
  }

  size_t fetchOtherId(size_t id) const
  {
    assert(popIndices_.size() == 2);

    if(popIndices_[0] == id)
      return popIndices_[1];

    else
      return popIndices_[0];
  }

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

  void cleanBasis_();

  // assigns two HetMoment pointers to each Pi2Moment (left and right loci)
  void linkPi2HetStats_();

  // exploits symmetry among statistics to reduce size of basis, given constraints imposed by selection
  void aliasMoments_();

  void compressBasis_();

};

#endif

