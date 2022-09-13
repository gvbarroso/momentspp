/*
 * Authors: Gustavo V. Barroso
 * Created: 05/08/2022
 * Last modified: 13/09/2022
 *
 */

#include "SumStatsLibrary.hpp"


void SumStatsLibrary::initStatsVector(size_t order = 2, const std::map<size_t, std::pair<size_t, size_t>>& popMap) // for a given epoch
{
  // map keys are populantion indices (vals are parents in previous epochs)
  order_ = 2;
  numPops_ = popMap.size();

  includeHetStats_(popMap);
  includeLdStats_(popMap);
}

void SumStatsLibrary::includeHetStats_(const std::map<size_t, std::pair<size_t, size_t>>& popMap)
{
  for(auto itI = std::begin(popMap); itI != std::end(popMap); ++itI)
    for(auto itJ = std::begin(popMap); itJ != std::end(popMap); ++itJ)
      stats_.try_emplace("H_" + bpp::TextTools::toString(itI->first) + bpp::TextTools::toString(itJ->first), 0.);
}

void SumStatsLibrary::includeLdStats_(const std::map<size_t, std::pair<size_t, size_t>>& popMap)
{
  for(auto itI = std::begin(popMap); itI != std::end(popMap); ++itI)
    for(auto itJ = std::begin(popMap); itJ != std::end(popMap); ++itJ)
      stats_.try_emplace("DD_" + bpp::TextTools::toString(itI->first) + bpp::TextTools::toString(itJ->first), 0.);

  for(auto itI = std::begin(popMap); itI != std::end(popMap); ++itI)
    for(auto itJ = std::begin(popMap); itJ != std::end(popMap); ++itJ)
      for(auto itK = std::begin(popMap); itK != std::end(popMap); ++itK)
        stats_.try_emplace("Dz_" + bpp::TextTools::toString(itI->first) + bpp::TextTools::toString(itJ->first) + bpp::TextTools::toString(itK->first), 0.);

  for(auto itI = std::begin(popMap); itI != std::end(popMap); ++itI)
    for(auto itJ = std::begin(popMap); itJ != std::end(popMap); ++itJ)
      for(auto itK = std::begin(popMap); itK != std::end(popMap); ++itK)
        for(auto itL = std::begin(popMap); itL != std::end(popMap); ++itL)
          stats_.try_emplace("pi2_" + bpp::TextTools::toString(itI->first) + bpp::TextTools::toString(itJ->first) + ";" + bpp::TextTools::toString(itK->first) + bpp::TextTools::toString(itL->first), 0.);
}
