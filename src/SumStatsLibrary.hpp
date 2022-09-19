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


// intent is to have one instance of SumStatsLibrary per Epoch
class SumStatsLibrary
{

private:
  size_t order_; // of moments
  size_t numPops_;

  std::vector<size_t> popIndices_; // independent object from popsMap_ for faster bookkeeping
  std::map<size_t, std::pair<size_t, size_t>> popsMap_; // pop-index -> pair of parents in prev epoch

  // WARNING split string on key doesn't work as intended when there are 2-digit population indices
  std::map<std::string, double> stats_; // stat name -> value NOTE use HashTable here?
  //https://dev.to/muiz6/string-hashing-in-c-1np3
  //https://stackoverflow.com/questions/8029121/how-to-hash-stdstring

public:
  SumStatsLibrary(size_t order = 2, const std::map<size_t, std::pair<size_t, size_t>>& popMap):
  order_(order),
  numPops_(popMap.size()),
  numDDStats_(0),
  numDDStats_(0),
  numDDStats_(0),
  numDDStats_(0),
  popIndices_(0),
  popsMap_(popMap),
  stats_(),
  popStatsMap_()
  {
    initStatsVector_();
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

  size_t getNumStats() const
  {
    return stats_.size();
  }

  const std::map<size_t, std::pair<size_t, size_t>>& getPopsMap() const
  {
    return popsMap_;
  }

  const std::map<std::string, double>& getStats() const
  {
    return stats_;
  }

  void setStatValue(const std::string& name, double value)
  {
    stats_.at(name) = value;
  }

  Eigen::VectorXd fetchYvec();

  bool hasPopIndex(const std::string& target, const std::string& query)
  {
    return countInstances(target, query) > 0;
  }

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
  void selectPopIndices_();

  void initStatsVector_();

  void includeHetStats_(const std::map<size_t, std::pair<size_t, size_t>>& popMap);

  void includeLdStats_(const std::map<size_t, std::pair<size_t, size_t>>& popMap);

  void compress_(); // exploits symmetry among statistics to reduce dimension of stats_

};

#endif

