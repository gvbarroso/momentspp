/*
 * Authors: Gustavo V. Barroso
 * Created: 05/08/2022
 * Last modified: 16/08/2022
 *
 */

#include "SumStatsLibrary.hpp"

#define  asString

void SumStatsLibrary::init(size_t numPops)
{
  setNumPops(numPops);

  includeHetStats_();
  includeLdStats_();
}

void SumStatsLibrary::includeHetStats_()
{
  for(size_t i = 0; i < numPops_; ++i)
    for(size_t j = 0; j < numPops_; ++j)
      stats_.try_emplace("H_" + asString_(i) + asString_(j), 0.);
}

void SumStatsLibrary::includeLdStats_()
{
  for(size_t i = 0; i < numPops_; ++i)
    for(size_t j = 0; j < numPops_; ++j)
      stats_.try_emplace("DD_" + asString_(i) + asString_(j), 0.);

  for(size_t i = 0; i < numPops_; ++i)
    for(size_t j = 0; j < numPops_; ++j)
      for(size_t k = 0; k < numPops_; ++k)
        stats_.try_emplace("Dz_" +  asString_(i) + asString_(j) + asString_(k), 0.);

  for(size_t i = 0; i < numPops_; ++i)
    for(size_t j = 0; j < numPops_; ++j)
      for(size_t k = 0; k < numPops_; ++k)
        for(size_t l = k; l < numPops_; ++l)
          if((k == i == l && j != i) && l < j)
            stats_.try_emplace("pi2_" + asString_(i) + asString_(j) + ";" + asString_(k) + asString_(l));
}
