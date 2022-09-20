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

void SumStatsLibrary::initMoments_()
{
  for(auto itI = std::begin(popIndices_); itI != std::end(popIndices_); ++itI)
  {
    for(auto itJ = std::begin(popIndices_); itJ != std::end(popIndices_); ++itJ)
    {
      moments_.push_back(Moment("DD_" + asString(*itI) + asString(*itJ), 0.));
      moments_.push_back(Moment("H_" + asString(*itI) + asString(*itJ), 0.));

      for(auto itK = std::begin(popIndices_); itK != std::end(popIndices_); ++itK)
      {
        moments_.push_back(Moment("Dz_" + asString(*itI) + asString(*itJ) + asString(*itK), 0.));

        for(auto itL = std::begin(popIndices_); itL != std::end(popIndices_); ++itL)
            moments_.push_back(Moment("pi2_" + asString(*itI) + asString(*itJ) + asString(*itK) + asString(*itL), 0.));
      }
    }
  }

  // NOTE this crucially determines the (lexicographical) order of stats in rows of Y and rows/cols of transition matrices
  std::sort(std::begin(moments_), std::end(moments_), [&](int a, int b) { a.getName() > b.getName() });
}
