/*
 * Authors: Gustavo V. Barroso
 * Created: 05/08/2022
 * Last modified: 19/09/2022
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

  std::vector<Moment> moments_; // independent object from popsMap_ for faster bookkeeping

public:
  SumStatsLibrary(size_t order = 2, const std::vector<size_t>& popIndices):
  order_(order),
  numPops_(popIndices.size()),
  moments_(0)
  {
    initMoments_(popIndices);
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

  void setHetMomentValue(size_t id1, size_t id2, double value)
  {
    std::string prefix = "H"; // just as a reminder
    // TODO set this efficiently by jumping over moments_
  }

  void setDdMomentValue(size_t id1, size_t id2, double value)
  {
    std::string prefix = "DD"; // just as a reminder
    // TODO set this efficiently by jumping over moments_
  }

  void setDzMomentValue(size_t id1, size_t id2, size_t id3, double value)
  {
    std::string prefix = "Dz"; // just as a reminder
    // TODO set this efficiently by jumping over moments_
  }

  void setPi2MomentValue(size_t id1, size_t id2, size_t id3, size_t id4, double value)
  {
    std::string prefix = "pi2"; // just as a reminder
    // TODO set this efficiently by jumping over moments_
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

  size_t indexLookup(const std::string& moment) const
  {
    return std::distance(std::begin(stats_), stats_.find(moment));
  }

private:
  void initMoments_(const std::vector<size_t>& popIndices);

  void includeHetStats_(const std::map<size_t, std::pair<size_t, size_t>>& popMap);

  void includeLdStats_(const std::map<size_t, std::pair<size_t, size_t>>& popMap);

  void compress_(); // exploits symmetry among statistics to reduce dimension of stats_

};

#endif

