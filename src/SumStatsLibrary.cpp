/*
 * Authors: Gustavo V. Barroso
 * Created: 05/08/2022
 * Last modified: 13/12/2022
 *
 */

#include <ios>

#include "SumStatsLibrary.hpp"

Eigen::VectorXd SumStatsLibrary::fetchYvec()
{
  Eigen::VectorXd y(moments_.size());

  for(size_t i = 0; i < moments_.size(); ++i)
    y(i) = moments_[i].getValue();

  return y;
}

void SumStatsLibrary::printMoments(std::ostream& stream)
{
  for(size_t i = 0; i < moments_.size(); ++i)
    stream << std::scientific << moments_[i].getPosition() << " | " << moments_[i].getName() << " = " << moments_[i].getValue() << "\n";
}

void SumStatsLibrary::initMoments_()
{
  for(auto itI = std::begin(popIndices_); itI != std::end(popIndices_); ++itI)
  {
    for(auto itJ = std::begin(popIndices_); itJ != std::end(popIndices_); ++itJ)
    {
      moments_.push_back(Moment("DD_" + asString(*itI) + "_" + asString(*itJ), 0.));

      // introduce ordered heterozigosity statistics Hp_ij (when i > j), Hq_ij (when i < j), and Hp_ii as well as Hq_ii
      if(*itI == *itJ) // if within population heterozigosity
      {
        moments_.push_back(Moment("Hp_" + asString(*itI) + "_" + asString(*itJ), 0.)); // p(1-p)
        moments_.push_back(Moment("Hq_" + asString(*itI) + "_" + asString(*itJ), 0.)); // (1-p)p
      }

      else if(*itI > *itJ) // first pop index higher
        moments_.push_back(Moment("Hp_" + asString(*itI) + "_" + asString(*itJ), 0.)); // p(1-p)

      else if(*itI < *itJ)
        moments_.push_back(Moment("Hq_" + asString(*itI) + "_" + asString(*itJ), 0.)); // (1-p)p

      for(auto itK = std::begin(popIndices_); itK != std::end(popIndices_); ++itK)
      {
        moments_.push_back(Moment("Dz_" + asString(*itI) + "_" + asString(*itJ) + "_" + asString(*itK), 0.));

        for(auto itL = std::begin(popIndices_); itL != std::end(popIndices_); ++itL)
          moments_.push_back(Moment("pi2_" + asString(*itI) + "_" + asString(*itJ) + "_" + asString(*itK) + "_" + asString(*itL), 0.));
      }
    }
  }

  // adds Dummy Moment lexicographically after H_ stats to convert into a homogeneous system (see Mutation::setUpMatrices_())
  moments_.push_back(Moment("I", 1.));

  // this determines the ascending lexicographical order of stats in the rows of transition matrices
  std::sort(std::begin(moments_), std::end(moments_), [](const Moment& a, const Moment& b) { return a.getName() < b.getName(); } );

  for(size_t i = 0; i < moments_.size(); ++i)
    moments_[i].setPosition(i);
}
