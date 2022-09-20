/*
 * Authors: Gustavo V. Barroso
 * Created: 05/08/2022
 * Last modified: 20/09/2022
 *
 */


#ifndef _SUMSTATSLIBRARY_H_
#define _SUMSTATSLIBRARY_H_

#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <cstring>
#include <utility>

#include <boost/algorithm/string.hpp>

#include <Bpp/Text/TextTools.h>

# include "Moment.hpp"

// intent is to have one instance of SumStatsLibrary per Epoch
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
  std::vector<Moment> moments_; // NOTE sorted lexicographically based on name_

public:
  SumStatsLibrary(size_t order = 2, const std::vector<size_t>& popIndices):
  order_(order),
  numPops_(popIndices.size()),
  numDDStats_(numPops_ * numPops_),
  numDzStats_(numPops_ * numPops_ * numPops_),
  numHetStats_(numPops_ * numPops_),
  numPi2Stats_(numPops_ * numPops_ * numPops_ * numPops_),
  popIndices_(popIndices),
  moments_(0)
  {
    std::sort(std::begin(popIndices_), std::end(popIndices_)); // safety
    initMoments_();
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

  const std::vector<Moment>& getMoments() const
  {
    return moments_;
  }

  size_t getNumStats() const
  {
    return moments_.size();
  }

  void setDdMomentValue(size_t id1, size_t id2, double value)
  {
    size_t rank1 = findPopIndexRank(id1);
    size_t rank2 = findPopIndexRank(id2);

    size_t focalMomIndex = findDdIndex(rank1, rank2);

    moments_[focalMomIndex].setValue(value);
  }

  void setDzMomentValue(size_t id1, size_t id2, size_t id3, double value)
  {
    size_t rank1 = findPopIndexRank(id1);
    size_t rank2 = findPopIndexRank(id2);
    size_t rank3 = findPopIndexRank(id3);

    size_t focalMomIndex = findDzIndex(rank1, rank2, rank3);

    moments_[focalMomIndex].setValue(value);
  }

  void setHetMomentValue(size_t id1, size_t id2, double value)
  {
    size_t rank1 = findPopIndexRank(id1);
    size_t rank2 = findPopIndexRank(id2);

    size_t focalMomIndex = findHetIndex(rank1, rank2);

    moments_[focalMomIndex].setValue(value);
  }

  void setPi2MomentValue(size_t id1, size_t id2, size_t id3, size_t id4, double value)
  {
     size_t rank1 = findPopIndexRank(id1);
     size_t rank2 = findPopIndexRank(id2);
     size_t rank3 = findPopIndexRank(id3);
     size_t rank4 = findPopIndexRank(id4);

     size_t focalMomIndex = findPi2Index(rank1, rank2, rank3, rank4);

     moments_[focalMomIndex].setValue(value);
  }

  size_t findPopIndexRank(size_t index) // among all pop indices
  {
    return std::distance(std::begin(popIndices_), std::binary_search(std::begin(popIndices_), std::end(popIndices_), index));
  }

  size_t findDdIndex(size_t rank1, size_t rank2)
  {
    return rank1 * numPops_ + rank2;
  }

  size_t findDzIndex(size_t rank1, size_t rank2, size_t rank3)
  {
    return numDDStats_ + rank1 * numPops_ * numPops_ + rank2 * numPops_ + rank3;
  }

  size_t findHetIndex(size_t rank1, size_t rank2)
  {
    return numDDStats_ + numDzStats_ + rank1 * numPops_ + rank2;
  }

  size_t findPi2Index(size_t rank1, size_t rank2, size_t rank3, size_t rank4)
  {
    return numDDStats_ + numDzStats_ + numHetStats_ + rank1 * numPops_ * numPops_ * numPops_ + rank2 * numPops_ * numPops_ + rank3 * numPops_ + rank4;
  }

  std::string asString(size_t i)
  {
    return bpp::TextTools::toString(i);
  }

  Eigen::VectorXd fetchYvec();

  std::vector<std::string> splitString(const std::string& target, const std::string& query) const
  {
    std::vector<std::string> ret(0);
    boost::split(ret, target, boost::is_any_of(query));

    return ret;
  }

  size_t countInstances(const std::string& target, const std::string& query) const
  {
    std::string::difference_type count = std::count(std::begin(target), std::end(target), *query.c_str());
    return static_cast<size_t>(count);
  }

  size_t indexLookup(const std::string& name) const
  {
    return std::distance(std::begin(moments_), moments_.find(name));
  }

private:
  void initMoments_();

  void compress_(); // exploits symmetry among statistics to reduce dimension of stats_

};

#endif

