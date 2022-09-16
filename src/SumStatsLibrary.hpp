/*
 * Authors: Gustavo V. Barroso
 * Created: 05/08/2022
 * Last modified: 16/09/2022
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
  size_t numDDStats_;
  size_t numDzStats_;
  size_t numHetStats_;
  size_t numPi2Stats_;

  std::vector<size_t> popIndices_; // independent object from popsMap_ for faster bookkeeping
  std::map<size_t, std::pair<size_t, size_t>> popsMap_; // pop-index -> pair of parents in prev epoch
  std::map<std::string, double> stats_; // stat name -> value

  // pop-index -> 4 vectors, each with indices of here pop-index appears in DD, Dz, Het and Pi2 stats, respectively
  // this object exists for fast access of this information when copying summary statistics (see Model::computeExpectedSumStats_()
  std::map<size_t, std::array<std::vector<size_t>, 4>> popStatsMap_;

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
    numDDStats_ = numPops_ * numPops_;
    numDzStats_ = numPops_ * numPops_ * numPops_;
    numHetStats_ = numPops_ * numPops_;
    numPi2Stats_ = numPops_ * numPops_ * numPops_ * numPops_;

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

  std::map<std::string, double>& getStats()
  {
    return stats_;
  }

  Eigen::VectorXd fetchYvec();

  void copyStatsToMap(Eigen::VectorXd& y);

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

  size_t findStatIndex(const std::string& mom); // fast version


private:
  void selectPopIndices_();

  void initStatsVector_();

  void includeHetStats_(const std::map<size_t, std::pair<size_t, size_t>>& popMap);

  void includeLdStats_(const std::map<size_t, std::pair<size_t, size_t>>& popMap);

  void compress_(); // exploits symmetry among statistics to reduce dimension of stats_

};

#endif

