/*
 * Authors: Gustavo V. Barroso
 * Created: 05/08/2022
 * Last modified: 15/09/2022
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
  bool compressed_; // indicates whether symmetrical statistics have been pulled together
  size_t order_; // of moments

  std::map<size_t, std::pair<size_t, size_t>> popsMap_; // pop index -> pair of parents in prev epoch
  std::map<std::string, double> stats_; // stat name -> value

public:
  SumStatsLibrary(size_t order = 2, const std::map<size_t, std::pair<size_t, size_t>>& popMap):
  compressed_(false),
  order_(order),
  popsMap_(popMap),
  stats_()
  {
    initStatsVector_();
  }

public:
  bool compressed()
  {
    return compressed_;
  }

  size_t getNumPops()
  {
    return popsMap_.size();
  }

  size_t getOrder()
  {
    return order_;
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

  std::map<std::string, double>& getStats()
  {
    return stats_;
  }

  std::vector<size_t> fetchPopIndices();

  void copyStats(Eigen::VectorXd& y);

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
  void initStatsVector_();

  void includeHetStats_(const std::map<size_t, std::pair<size_t, size_t>>& popMap);

  void includeLdStats_(const std::map<size_t, std::pair<size_t, size_t>>& popMap);

  void compress_(); // exploits symmetry among statistics to reduce dimension of stats_

};

#endif

