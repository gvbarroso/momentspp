/*
 * Authors: Gustavo V. Barroso
 * Created: 05/08/2022
 * Last modified: 09/02/2023
 *
 */


#ifndef _SUMSTATSLIBRARY_H_
#define _SUMSTATSLIBRARY_H_

#include <string>
#include <vector>
#include <memory>
#include <algorithm>
#include <cstring>
#include <utility>
#include <ostream>

#include <Eigen/Core>

#include <boost/algorithm/string.hpp>

#include <Bpp/Text/TextTools.h>

#include "Population.hpp"
#include "Moment.hpp"
#include "DdMoment.hpp"
#include "DzMoment.hpp"
#include "HetMoment.hpp"
#include "Pi2Moment.hpp"


// intent is to have one instance of SumStatsLibrary per Epoch because each Epoch potentially has unique population sets
class SumStatsLibrary
{

private:
  size_t order_; // of moments
  size_t numPops_;
  size_t numDDStats_;
  size_t numDzStats_;
  size_t numHetStats_;
  size_t numPi2Stats_;

  std::vector<size_t> popIndices_; // among all Moments, stored for bookkeeping
  std::vector<std::string> hetSuffixes_;
  std::vector<std::string> pi2Suffixes_;
  std::vector<std::shared_ptr<Moment>> moments_; // sorted lexicographically based on their name_

public:
  SumStatsLibrary():
  order_(0),
  numPops_(0),
  numDDStats_(0),
  numDzStats_(0),
  numHetStats_(0),
  numPi2Stats_(0),
  popIndices_(0),
  hetSuffixes_(0),
  pi2Suffixes_(0),
  moments_(0)
  { }

  SumStatsLibrary(size_t order, const std::map<size_t, std::shared_ptr<Population>>& popMap):
  order_(order),
  numPops_(popMap.size()),
  numDDStats_(numPops_ * numPops_),
  numDzStats_(numPops_ * numPops_ * numPops_),
  numHetStats_(2 * numPops_ * numPops_), // inits with all sampling permutations p(1-p), (1-p)p
  numPi2Stats_(4 * numPops_ * numPops_ * numPops_ * numPops_), // inits with all sampling permutations p(1-p), (1-p)p for each locus
  popIndices_(0),
  hetSuffixes_(0),
  pi2Suffixes_(0),
  moments_(0)
  {
    hetSuffixes_ = { "A", "B" }; // p(1-p), (1-p)p
    pi2Suffixes_ = { "A", "B", "C", "D" }; // p(1-p)*(1-q)q, p(1-p)*q(1-q), (1-p)1*(1-q)q, (1-p)p*q(1-q)

    popIndices_.reserve(popMap.size());

    for(auto it = std::begin(popMap); it != std::end(popMap); ++it)
    {
      assert(it->first == it->second->getId());
      popIndices_.emplace_back(it->first);
    }

    std::sort(std::begin(popIndices_), std::end(popIndices_));
    initMoments_(popMap);
  }

public:
  size_t getNumPops()
  {
    return numPops_;
  }

  size_t getOrder()
  {
    return order_;
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

  size_t getNumDDStats() const
  {
    return numDDStats_;
  }

  size_t getNumDzStats() const
  {
    return numDzStats_;
  }

  size_t getNumHetStats() const
  {
    return numHetStats_;
  }

  size_t getNumPi2Stats() const
  {
    return numPi2Stats_;
  }

  size_t getNumStats() const
  {
    return 1 + numDDStats_ + numDzStats_ + numHetStats_ + numPi2Stats_;
  }

  // naive search
  std::shared_ptr<Moment> getMoment(const std::string& name) const
  {
    std::shared_ptr<Moment> ptr = nullptr; // object slicing, can only get (base) Moment attributes
    for(auto itMom = std::begin(moments_); itMom != std::end(moments_); ++itMom)
    {
      if((*itMom)->getName() == name)
        ptr = *itMom;
    }

    assert(ptr != nullptr);
    return ptr;
  }

  // the following methods use pop-ids to track down moments' positions inside moments_ vector (see Model::linkMoments_())

  std::shared_ptr<DdMoment> getDdMoment(size_t id1, size_t id2) const
  {
    size_t focalMomIndex = findDdIndex(id1, id2);
    auto ret = std::dynamic_pointer_cast<DdMoment>(moments_[focalMomIndex]);

    if(ret != nullptr)
      return ret;

    else
      throw bpp::Exception("SumStatsLibrary::bad dynamic_pointer_cast attempt!");
  }

  std::shared_ptr<DzMoment> getDzMoment(size_t id1, size_t id2, size_t id3) const
  {
    size_t focalMomIndex = findDzIndex(id1, id2, id3);
    auto ret = std::dynamic_pointer_cast<DzMoment>(moments_[focalMomIndex]);

    if(ret != nullptr)
      return ret;

    else
      throw bpp::Exception("SumStatsLibrary::bad dynamic_pointer_cast attempt!");
  }

  std::shared_ptr<HetMoment> getHetMoment(size_t id1, size_t id2, const std::string& suffix) const
  {
    size_t focalMomIndex = findHetIndex(id1, id2, suffix);
    auto ret = std::dynamic_pointer_cast<HetMoment>(moments_[focalMomIndex]);

    if(ret != nullptr)
      return ret;

    else
      throw bpp::Exception("SumStatsLibrary::bad dynamic_pointer_cast attempt!");
  }

  std::shared_ptr<Pi2Moment> getPi2Moment(size_t id1, size_t id2, size_t id3, size_t id4, const std::string& suffix) const
  {
    size_t focalMomIndex = findPi2Index(id1, id2, id3, id4, suffix);
    auto ret = std::dynamic_pointer_cast<Pi2Moment>(moments_[focalMomIndex]);

    if(ret != nullptr)
      return ret;

    else
      throw bpp::Exception("SumStatsLibrary::bad dynamic_pointer_cast attempt!");
  }

  std::shared_ptr<Moment> getDummyMoment() const
  {
    return moments_[getDummyIndex()];
  }

  size_t findPopIndexRank(size_t index) const // among all pop indices
  {
    auto it = std::lower_bound(std::begin(popIndices_), std::end(popIndices_), index);
    assert(it != std::end(popIndices_));

    return std::distance(std::begin(popIndices_), it); // indexed from 0
  }

  size_t findHetSuffixRank(const std::string& suffix) const
  {
    auto it = std::lower_bound(std::begin(hetSuffixes_), std::end(hetSuffixes_), suffix);
    assert(it != std::end(hetSuffixes_));

    return std::distance(std::begin(hetSuffixes_), it); // indexed from 0
  }

  size_t findPi2SuffixRank(const std::string& suffix) const
  {
    auto it = std::lower_bound(std::begin(pi2Suffixes_), std::end(pi2Suffixes_), suffix);
    assert(it != std::end(pi2Suffixes_));

    return std::distance(std::begin(pi2Suffixes_), it); // indexed from 0
  }

  size_t findDdIndex(size_t id1, size_t id2) const
  {
    size_t rank1 = findPopIndexRank(id1);
    size_t rank2 = findPopIndexRank(id2);

    return rank1 * numPops_ + rank2;
  }

  size_t findDzIndex(size_t id1, size_t id2, size_t id3) const
  {
    size_t rank1 = findPopIndexRank(id1);
    size_t rank2 = findPopIndexRank(id2);
    size_t rank3 = findPopIndexRank(id3);

    return numDDStats_ + rank1 * numPops_ * numPops_ + rank2 * numPops_ + rank3;
  }

  size_t findHetIndex(size_t id1, size_t id2, const std::string& suffix) const
  {
    size_t rank1 = findPopIndexRank(id1);
    size_t rank2 = findPopIndexRank(id2);

    size_t suffixRank = findHetSuffixRank(suffix);

    return numDDStats_ + numDzStats_ + rank1 * numPops_ * hetSuffixes_.size() + rank2 * hetSuffixes_.size() + suffixRank;
  }

  size_t getDummyIndex() const
  {
    return numDDStats_ + numDzStats_ + numHetStats_;
  }

  size_t findPi2Index(size_t id1, size_t id2, size_t id3, size_t id4, const std::string& suffix) const
  {
    size_t rank1 = findPopIndexRank(id1);
    size_t rank2 = findPopIndexRank(id2);
    size_t rank3 = findPopIndexRank(id3);
    size_t rank4 = findPopIndexRank(id4);

    size_t suffixRank = findPi2SuffixRank(suffix);

    // 1 + because of dummy Moment "I_" after "H_**" to make system homogeneous(see initMoments_())
    return 1 + numDDStats_ + numDzStats_ + numHetStats_ + rank1 * numPops_ * numPops_ * numPops_ * pi2Suffixes_.size() + rank2 * numPops_ * numPops_ * pi2Suffixes_.size() + rank3 * numPops_ * pi2Suffixes_.size() + rank4 * pi2Suffixes_.size() + suffixRank;
  }

  std::string asString(size_t i)
  {
    return bpp::TextTools::toString(i);
  }

  Eigen::VectorXd fetchYvec();

  void printMoments(std::ostream& stream);

  // exploits symmetry among statistics to reduce dimension of stats_, given constraints imposed by Selection
  void aliasMoments(const std::vector<size_t>& selectedPopIds);

private:
  void initMoments_(const std::map<size_t, std::shared_ptr<Population>>& popMap);

  // assign two HetMoment pointers to each Pi2Moment (left and right loci)
  void linkPi2HetStats_();

};

#endif

