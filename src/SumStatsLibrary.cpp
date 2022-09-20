/*
 * Authors: Gustavo V. Barroso
 * Created: 05/08/2022
 * Last modified: 20/09/2022
 *
 */


#include "SumStatsLibrary.hpp"


Eigen::VectorXd SumStatsLibrary::fetchYvec()
{
  Eigen::VectorXd& y(moments_.size());

  for(size_t i = 0; i < moments_.size(); ++i)
    y(0, i) = moments_[i].getValue();

  return y;
}

// NOTE this method crucially determines the order of stats in rows of Y and rows/cols of transition matrices
void SumStatsLibrary::initMoments_()
{
  includeHetStats_(popIndices_);
  includeLdStats_(popIndices_);
}

void SumStatsLibrary::includeHetStats_()
{
  for(auto itI = std::begin(popIndices_); itI != std::end(popIndices_); ++itI)
    for(auto itJ = std::begin(popIndices_); itJ != std::end(popIndices_); ++itJ)
      moments_.push_back(Moment("H_" + bpp::TextTools::toString(*itI) + bpp::TextTools::toString(*itJ), 0.));
}

void SumStatsLibrary::includeLdStats_()
{
  for(auto itI = std::begin(popIndices_); itI != std::end(popIndices_); ++itI)
    for(auto itJ = std::begin(popMap); itJ != std::end(popIndices_); ++itJ)
      moments_.push_back(Moment("DD_" + bpp::TextTools::toString(*itI) + bpp::TextTools::toString(*itJ), 0.));

  for(auto itI = std::begin(popIndices_); itI != std::end(popIndices_); ++itI)
    for(auto itJ = std::begin(popIndices_); itJ != std::end(popIndices_); ++itJ)
      for(auto itK = std::begin(popIndices_); itK != std::end(popIndices_); ++itK)
        moments_.push_back(Moment("Dz_" + bpp::TextTools::toString(*itI) + bpp::TextTools::toString(*itJ) + bpp::TextTools::toString(*itK), 0.));

  for(auto itI = std::begin(popIndices_); itI != std::end(popIndices_); ++itI)
    for(auto itJ = std::begin(popIndices_); itJ != std::end(popIndices_); ++itJ)
      for(auto itK = std::begin(popIndices_); itK != std::end(popIndices_); ++itK)
        for(auto itL = std::begin(popIndices_); itL != std::end(popIndices_); ++itL)
          moments_.push_back(Moment("pi2_" + bpp::TextTools::toString(*itI) + bpp::TextTools::toString(*itJ) + ";" + bpp::TextTools::toString(*itK) + bpp::TextTools::toString(*itL), 0.));
}
