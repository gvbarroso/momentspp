/*
 * Authors: Gustavo V. Barroso
 * Created: 05/08/2022
 * Last modified: 05/08/2022
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
  std::map<std::string, double> stats_;

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
