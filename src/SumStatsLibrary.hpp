/*
 * Authors: Gustavo V. Barroso
 * Created: 05/08/2022
 * Last modified: 12/12/2022
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
#include <ostream>

#include <Eigen/Core>

#include <boost/algorithm/string.hpp>

#include <Bpp/Text/TextTools.h>

#include "Moment.hpp"
#include "Population.hpp"


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
  std::vector<Moment> moments_; // sorted lexicographically based on name_

public:
  SumStatsLibrary():
  order_(0),
  numPops_(0),
  numDDStats_(0),
  numDzStats_(0),
  numHetStats_(0),
  numPi2Stats_(0),
  popIndices_(0),
  moments_(0)
  { }

  SumStatsLibrary(size_t order, const std::map<size_t, std::shared_ptr<Population>>& popMap):
  order_(order),
  numPops_(popMap.size()),
  numDDStats_(numPops_ * numPops_),
  numDzStats_(numPops_ * numPops_ * numPops_),
  numHetStats_(numPops_ * numPops_ + numPops_ * (numPops_ - 1)), // sorted H statistics (including within pops)
  numPi2Stats_(numPops_ * numPops_ * numPops_ * numPops_),
  popIndices_(0),
  moments_(0)
  {
    popIndices_.reserve(popMap.size());

    for(auto it = std::begin(popMap); it != std::end(popMap); ++it)
      popIndices_.emplace_back(it->first);

    std::sort(std::begin(popIndices_), std::end(popIndices_));
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

  const std::vector<size_t>& getPopIndices() const
  {
    return popIndices_;
  }

  const std::vector<Moment>& getMoments() const
  {
    return moments_;
  }

  std::vector<Moment>& getMoments()
  {
    return moments_;
  }

  size_t getNumStats() const
  {
    return moments_.size();
  }

  // these methods use pop-ids to track down moments' positions inside moments_ vector (see Model::linkMoments_())

  const Moment& getDdMoment(size_t id1, size_t id2) const
  {
    size_t focalMomIndex = findDdIndex(id1, id2);
    return moments_[focalMomIndex];
  }

  Moment& getDdMoment(size_t id1, size_t id2)
  {
    size_t focalMomIndex = findDdIndex(id1, id2);
    return moments_[focalMomIndex];
  }

  const Moment& getDzMoment(size_t id1, size_t id2, size_t id3) const
  {
    size_t focalMomIndex = findDzIndex(id1, id2, id3);
    return moments_[focalMomIndex];
  }

  Moment& getDzMoment(size_t id1, size_t id2, size_t id3)
  {
    size_t focalMomIndex = findDzIndex(id1, id2, id3);
    return moments_[focalMomIndex];
  }

  const Moment& getHetMoment(size_t id1, size_t id2) const
  {
    size_t focalMomIndex = findHetIndex(id1, id2);
    return moments_[focalMomIndex];
  }

  Moment& getHetMoment(size_t id1, size_t id2)
  {
    size_t focalMomIndex = findHetIndex(id1, id2);
    return moments_[focalMomIndex];
  }

  const Moment& getPi2Moment(size_t id1, size_t id2, size_t id3, size_t id4) const
  {
    size_t focalMomIndex = findPi2Index(id1, id2, id3, id4);
    return moments_[focalMomIndex];
  }

  Moment& getPi2Moment(size_t id1, size_t id2, size_t id3, size_t id4)
  {
    size_t focalMomIndex = findPi2Index(id1, id2, id3, id4);
    return moments_[focalMomIndex];
  }

  const Moment& getDummyMoment() const
  {
    return moments_[getDummyIndex()];
  }

  Moment& getDummyMoment()
  {
    return moments_[getDummyIndex()];
  }

  void setDdMomentValue(size_t id1, size_t id2, double value)
  {
    getDdMoment(id1, id2).setValue(value);
  }

  void setDzMomentValue(size_t id1, size_t id2, size_t id3, double value)
  {
    getDzMoment(id1, id2, id3).setValue(value);
  }

  void setHetMomentValue(size_t id1, size_t id2, double value)
  {
    getHetMoment(id1, id2).setValue(value);
  }

  void setPi2MomentValue(size_t id1, size_t id2, size_t id3, size_t id4, double value)
  {
    getPi2Moment(id1, id2, id3, id4).setValue(value);
  }

  size_t findPopIndexRank(size_t index) const // among all pop indices
  {
    return std::distance(std::begin(popIndices_), std::lower_bound(std::begin(popIndices_), std::end(popIndices_), index)); // indexed from 0
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

  size_t findHetIndex(size_t id1, size_t id2) const // WARNING
  {
    size_t rank1 = findPopIndexRank(id1);
    size_t rank2 = findPopIndexRank(id2);

    return numDDStats_ + numDzStats_ + rank1 * numPops_ + rank2 * numPops_; // TODO correct for ordered stats
  }

  size_t findPi2Index(size_t id1, size_t id2, size_t id3, size_t id4) const
  {
    size_t rank1 = findPopIndexRank(id1);
    size_t rank2 = findPopIndexRank(id2);
    size_t rank3 = findPopIndexRank(id3);
    size_t rank4 = findPopIndexRank(id4);

    // NOTE 1 + because of dummy Moment "I_" after "H_**" to make system homogeneous(see initMoments_())
    return 1 + numDDStats_ + numDzStats_ + numHetStats_ + rank1 * numPops_ * numPops_ * numPops_ + rank2 * numPops_ * numPops_ + rank3 * numPops_ + rank4;
  }

  size_t getDummyIndex() const
  {
    return numDDStats_ + numDzStats_ + numHetStats_;
  }

  std::string asString(size_t i)
  {
    return bpp::TextTools::toString(i);
  }

  Eigen::VectorXd fetchYvec();

  void printMoments(std::ostream& stream);

private:
  void initMoments_();

  void compress_(); // exploits symmetry among statistics to reduce dimension of stats_

};

#endif

