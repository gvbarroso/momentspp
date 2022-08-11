/*
 * Authors: Gustavo V. Barroso
 * Created: 05/08/2022
 * Last modified: 05/08/2022
 *
 */

#include "SumStatsLibrary.hpp"

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
      stats_.try_emplace("H_" + bpp::TextTools::toString(i) + "_" + bpp::TextTools::toString(j), 0.);
}

void SumStatsLibrary::includeLdStats_()
{
  for(size_t i = 0; i < numPops_; ++i)
    for(size_t j = 0; j < numPops_; ++j)
      stats_.try_emplace("DD_" + bpp::TextTools::toString(i) + "_" + bpp::TextTools::toString(j), 0.);

  for(size_t i = 0; i < numPops_; ++i)
    for(size_t j = 0; j < numPops_; ++j)
      for(size_t k = 0; k < numPops_; ++k)
        stats_.try_emplace("D_" +  bpp::TextTools::toString(i) +
                           "_z_" + bpp::TextTools::toString(j), "_",
                                   bpp::TextTools::toString(k), 0.);

  for(size_t i = 0; i < numPops_; ++i)
    for(size_t j = 0; j < numPops_; ++j)
      for(size_t k = 0; k < numPops_; ++k)
        for(size_t l = k; l < numPops_; ++l)
          if((k == i == l && j != i) && l < j)
            stats_.try_emplace("pi2_" + bpp::TextTools::toString(i) + "_" +
                                        bpp::TextTools::toString(j) + "_" +
                                        bpp::TextTools::toString(k) + "_" +
                                        bpp::TextTools::toString(l));
}
