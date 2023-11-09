/*
 * Authors: Gustavo V. Barroso
 * Created: 09/11/2022
 * Last modified: 09/11/2023
 */


#ifndef _SIM_TOOLS_H_
#define _SIM_TOOLS_H_

#include <set>
#include <string>
#include <vector>
#include <algorithm>
#include <cstring>
#include <utility>
#include <ostream>
#include <cassert>

#include "TwoLocusPop.hpp"

// main purpose is to compute cross-population summary statistics
class SimTools
{

public:
  SimTools():
  { }

public:
  double fetchD()
  {
    double n = getNumHaps();
    return (count_ab_ * count_AB_ - count_Ab_ * count_aB_) / (n * n);
  }

  double fetchHl()
  {
    return fetchP() * (1. - fetchP());
  }

  double fetchHr()
  {
    return fetchQ() * (1. - fetchQ());
  }

  double fetchDz()
  {
    return fetchD() * (1. - 2. * fetchP()) * (1. - 2. * fetchQ());
  }

  double fetchDsqr()
  {
    double n = getNumHaps();
    return ((count_ab_ * count_AB_ - count_Ab_ * count_aB_) * (count_ab_ * count_AB_ - count_Ab_ * count_aB_)) / (n * n * n * n);
  }

  double fetchPi2()
  {
    return fetchHl() * fetchHr();
  }

};

#endif

