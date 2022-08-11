/*
 * Authors: Gustavo V. Barroso
 * Created: 05/08/2022
 * Last modified: 09/08/2022
 *
 */


#ifndef _SUMSTATSLIBRARY_H_
#define _SUMSTATSLIBRARY_H_

#include <string>
#include <vector>
#include <map>

#include <Bpp/Text/TextTools.h>

#include "OptionsContainer.h"


class SumStatsLibrary
{

private:
  size_t numPops_;
  size_t order_; // of the moment; number of tracked lineages
  std::map<std::string, double> stats_;  // name -> value

public:
  SumStatsLibrary():
  numPops_(1)
  { }

  SumStatsLibrary(size_t numPops):
  numPops_(numPops)
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

private:
  void includeHetStats_();

  void includeLdStats_();

};
