/*
 * Authors: Gustavo V. Barroso
 * Created: 05/08/2022
 * Last modified: 12/08/2022
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
  std::map<std::string, double> stats_;  // name -> value

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

  size_t getNumPops()
  {
    return numPops;
  }

  void setNumPops(size_t numPops)
  {
    numPops_ = numPops;
  }

  size_t getOrder()
  {
    return order_;
  }

  size_t getNumStats()
  {
    return stats_.size()
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
    getNumDDStats() + getNumDzStats();
  }

  size_t getNumDDStats()
  {
    return numPops_ + (numPops_ * (numPops_ - 1) / 2); // P + (P choose 2) (works for order_ == 2)
  }

  size_t getNumDzStats()
  {
    return numPops_ * numPops_ * numPops_;
  }

  size_t getNumHetStats()
  {
    numPops_ + (numPops_ * (numPops_ - 1) / 2); // P + (P choose 2) (works for order_ == 2)
  }

  size_t getNumPiTwoStats()
  {
    (numPops_ + (numPops_ * (numPops_ - 1) / 2)) * (numPops_ + (numPops_ * (numPops_ - 1) / 2));
  }

  size_t getNumDiversityStats()
  {
    getNumHetStats() + getNumPiTwoStats();
  }

private:
  void includeHetStats_();

  void includeLdStats_();

  std::vector<std::string> splitUnder_(const std::string& s)
  {
    std::vector<std::string> ret;
    boost::split(ret, s, boost::is_any_of("_"));

    return ret;
  }

  std::string asString_(size_t i)
  {
    return bpp::TextTools::toString(i);
  }

  size_t countChars_(const std::string& s, char c)
  {
    std::string::difference_type n = std::count(std::begin(s), std::end(s), '_');
    return static_cast<size_t>(n);
  }

};



