/*
 * Authors: Gustavo V. Barroso
 * Created: 05/08/2022
 * Last modified: 14/09/2022
 *
 */

#include "SumStatsLibrary.hpp"


// NOTE this method crucially determines the order of stats in rows of Y and rows/cols of transition matrices
void SumStatsLibrary::initStatsVector(const std::map<size_t, std::pair<size_t, size_t>>& popMap) // for a given epoch
{
  // map keys are populantion indices (vals are parents in previous epochs)
  std::vector<size_t> popIndices = fetchPopIndices();

  includeHetStats_(popIndices);
  includeLdStats_(popIndices);
}

std::vector<size_t> SumStatsLibrary::fetchPopIndices()
{
  std::vector<size_t> popIndices(0);
  popIndices.reserve(popMap_.size());

  for(auto it = std::begin(popMap_); it != std::end(popMap_); ++it)
    popIndices.emplace_back(it->first);

  return popIndices;
}

void SumStatsLibrary::includeHetStats_(const std::vector<size_t>& popIndices)
{
  for(auto itI = std::begin(popIndices); itI != std::end(popIndices); ++itI)
    for(auto itJ = std::begin(popIndices); itJ != std::end(popIndices); ++itJ)
      stats_.try_emplace("H_" + bpp::TextTools::toString(*itI) + bpp::TextTools::toString(*itJ), 0.);
}

void SumStatsLibrary::includeLdStats_(const std::vector<size_t>& popIndices)
{
  for(auto itI = std::begin(popIndices); itI != std::end(popIndices); ++itI)
    for(auto itJ = std::begin(popMap); itJ != std::end(popIndices); ++itJ)
      stats_.try_emplace("DD_" + bpp::TextTools::toString(*itI) + bpp::TextTools::toString(*itJ), 0.);

  for(auto itI = std::begin(popIndices); itI != std::end(popIndices); ++itI)
    for(auto itJ = std::begin(popIndices); itJ != std::end(popIndices); ++itJ)
      for(auto itK = std::begin(popIndices); itK != std::end(popIndices); ++itK)
        stats_.try_emplace("Dz_" + bpp::TextTools::toString(*itI) + bpp::TextTools::toString(*itJ) + bpp::TextTools::toString(*itK), 0.);

  for(auto itI = std::begin(popIndices); itI != std::end(popIndices); ++itI)
    for(auto itJ = std::begin(popIndices); itJ != std::end(popIndices); ++itJ)
      for(auto itK = std::begin(popIndices); itK != std::end(popIndices); ++itK)
        for(auto itL = std::begin(popIndices); itL != std::end(popIndices); ++itL)
          stats_.try_emplace("pi2_" + bpp::TextTools::toString(*itI) + bpp::TextTools::toString(*itJ) + ";" + bpp::TextTools::toString(*itK) + bpp::TextTools::toString(*itL), 0.);
}
