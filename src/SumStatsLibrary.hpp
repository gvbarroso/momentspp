/*
 * Authors: Gustavo V. Barroso
 * Created: 05/08/2022
 * Last modified: 19/08/2022
 *
 */


#ifndef _SUMSTATSLIBRARY_H_
#define _SUMSTATSLIBRARY_H_

#include <string>
#include <vector>
#include <map>
#include <algorithm>

#include <boost/algorithm/string.hpp>

#include <Bpp/Text/TextTools.h>

#include "OptionsContainer.hpp"


class SumStatsLibrary
{

private:
  size_t numPops_;
  size_t order_; // of the moment; number of tracked lineages
  // TODO make this a std::map<std::string, std::pair<size_t, double>> stats_; where the size_t represents the redundancy factor of each statistic (>= 1)
  std::map<std::string, double> stats_;  // name -> value (Y vector)
  // NOTE row and column order of matrices follow lexicographical order of stats_'s names

public:
  SumStatsLibrary():
  numPops_(1),
  order_(2)
  { }

  SumStatsLibrary(size_t numPops):
  numPops_(numPops),
  order_(2)
  { }

  SumStatsLibrary(size_t order):
  numPops_(1),
  order_(order)
  { }

  SumStatsLibrary(size_t numPops, size_t order):
  numPops_(numPops),
  order_(order)
  { }

  SumStatsLibrary(const OptionsContainer& opt):
  numPops_(opt->getNumberOfPopulations()),
  order_(opt->getOrder())
  { }

  size_t getNumPops()
  {
    return numPops;
  }

  void setNumPops(const OptionsContainer& opt)
  {
    numPops_ = opt->getNumberOfPopulations();
  }

  size_t getOrder()
  {
    return order_;
  }

  void setOrder(const OptionsContainer& opt)
  {
    order_ = opt->getOrder();
  }

  size_t getNumStats()
  {
    return stats_.size();
  }

  const std::map<std::string, double>& getStats()
  {
    return stats_;
  }

  double getStat(const std::string& name)
  {
    return stats_.at(name);
  }

  void init();

  size_t getNumLdStats()
  {
    return getNumDDStats() + getNumDzStats();
  }

  size_t getNumDDStats()
  {
    return numPops_ + (numPops_ * (numPops_ - 1) / 2); // P + (P choose 2)
  }

  size_t getNumDzStats()
  {
    return numPops_ * numPops_ * numPops_;
  }

  size_t getNumHetStats()
  {
    return numPops_ + (numPops_ * (numPops_ - 1) / 2); // P + (P choose 2)
  }

  size_t getNumPiTwoStats()
  {
    return (numPops_ + (numPops_ * (numPops_ - 1) / 2)) * (numPops_ + (numPops_ * (numPops_ - 1) / 2));
  }

  size_t getNumDiversityStats()
  {
    return getNumHetStats() + getNumPiTwoStats();
  }

  std::vector<std::string> splitString(const std::string& target, const std::string& query)
  {
    std::vector<std::string> ret(0);
    boost::split(ret, target, boost::is_any_of(query));

    return ret;
  }

  std::string asString(size_t i)
  {
    return bpp::TextTools::toString(i);
  }

  size_t countInstances(const std::string& target, const std::string& query)
  {
    std::string::difference_type count = std::count(std::begin(target), std::end(target), query);
    return static_cast<size_t>(count);
  }

  size_t indexLookup(const std::string& moment)
  {
    auto pos = stats_.find(moment);

    if(pos == std::end(stats_)
      throw bpp::Exception("Moments++::SumStatsLibrary::Could not find index of moment " + moment);

    else
      return pos - std::begin(stats_);
  }

private:
  void includeHetStats_();

  void includeLdStats_();

  void compress_(); // exploits redundancy among statistics to reduce dimension

};



