/*
 * Authors: Gustavo V. Barroso
 * Created: 05/08/2022
 * Last modified: 16/09/2022
 *
 */


#include "SumStatsLibrary.hpp"


Eigen::VectorXd SumStatsLibrary::fetchYvec()
{
  Eigen::VectorXd& y(stats_.size());

  for(size_t i = 0; i < stats_.size(); ++i)
    y(0, i) = stats_[i].second;

  return y;
}

void SumStatsLibrary::copyStatsToMap(const Eigen::VectorXd& y)
{
  if(y.size() != stats_.size())
    throw bpp::Exception("SUM_STATS_LIBRARY::Attempted to copy from vector of different size!");

  else
    for(size_t i = 0; i < stats_.size(); ++i)
      stats_[i].second = y(0, i);
}

// TODO use the next 4 methods to fill in popStatsMap_
void fillDDStatsIndices_()
{
  size_t numPopIdDDStats = numDDstats_ - ((numPops_ - 1) * (numPops_ - 1)); // or just numPops + numPops - 1?


}

void fillDzStatsIndices_()
{
  size_t numPopIdDzStats = numDzStats_ - ((numPops_ - 1) * (numPops_ - 1) * (numPops_ - 1));

}

// returns indices of Het stats concerning index popId
void fillHetStatsIndices_()
{
  size_t numPopIdHetStats = numHetStats_ - ((numPops_ - 1) * (numPops_ - 1));

}

// returns indices of Pi2 stats concerning index popId
void fillPi2StatsIndices_()
{
  size_t numPopIdPi2Stats = numPi2Stats_ - ((numPops_ - 1) * (numPops_ - 1) * (numPops_ - 1) * (numPops_ - 1));

}

void SumStatsLibrary::selectPopIndices_()
{
  popIndices_.reserve(popsMap_.size());

  for(auto it = std::begin(popsMap_); it != std::end(popsMap_); ++it)
    popIndices_.emplace_back(it->first);
}

// NOTE this method crucially determines the order of stats in rows of Y and rows/cols of transition matrices
void SumStatsLibrary::initStatsVector_()
{
  // map keys are populantion indices (vals are parents in previous epochs)
  std::vector<size_t> popIndices = fetchPopIndices();

  includeHetStats_(popIndices);
  includeLdStats_(popIndices);
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
