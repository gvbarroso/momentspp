/*
 * Authors: Gustavo V. Barroso
 * Created: 05/08/2022
 * Last modified: 21/09/2022
 *
 */


#include "SumStatsLibrary.hpp"

const Moment& SumStatsLibrary::getDdMoment(size_t id1, size_t id2) const
{
  size_t rank1 = findPopIndexRank(id1);
  size_t rank2 = findPopIndexRank(id2);
  size_t focalMomIndex = findDdIndex(rank1, rank2);

  return moments_[focalMomIndex];
}

const Moment& SumStatsLibrary::getDzMoment(size_t id1, size_t id2, size_t id3) const
{
  size_t rank1 = findPopIndexRank(id1);
  size_t rank2 = findPopIndexRank(id2);
  size_t rank3 = findPopIndexRank(id3);
  size_t focalMomIndex = findDzIndex(rank1, rank2, rank3);

  return moments_[focalMomIndex];
}

const Moment& SumStatsLibrary::getHetMoment(size_t id1, size_t id2) const
{
  size_t rank1 = findPopIndexRank(id1);
  size_t rank2 = findPopIndexRank(id2);
  size_t focalMomIndex = findHetIndex(rank1, rank2);

  return moments_[focalMomIndex];
}

const Moment& SumStatsLibrary::getPi2Moment(size_t id1, size_t id2, size_t id3, size_t id4) const
{
  size_t rank1 = findPopIndexRank(id1);
  size_t rank2 = findPopIndexRank(id2);
  size_t rank3 = findPopIndexRank(id3);
  size_t rank4 = findPopIndexRank(id4);
  size_t focalMomIndex = findPi2Index(rank1, rank2, rank3, rank4);

  return moments_[focalMomIndex];
}

void SumStatsLibrary::setDdMomentValue(size_t id1, size_t id2, double value)
{
  size_t rank1 = findPopIndexRank(id1);
  size_t rank2 = findPopIndexRank(id2);
  size_t focalMomIndex = findDdIndex(rank1, rank2);

  moments_[focalMomIndex].setValue(value);
}

void SumStatsLibrary::setDzMomentValue(size_t id1, size_t id2, size_t id3, double value)
{
  size_t rank1 = findPopIndexRank(id1);
  size_t rank2 = findPopIndexRank(id2);
  size_t rank3 = findPopIndexRank(id3);
  size_t focalMomIndex = findDzIndex(rank1, rank2, rank3);

  moments_[focalMomIndex].setValue(value);
}

void SumStatsLibrary::setHetMomentValue(size_t id1, size_t id2, double value)
{
  size_t rank1 = findPopIndexRank(id1);
  size_t rank2 = findPopIndexRank(id2);
  size_t focalMomIndex = findHetIndex(rank1, rank2);

  moments_[focalMomIndex].setValue(value);
}

void SumStatsLibrary::setPi2MomentValue(size_t id1, size_t id2, size_t id3, size_t id4, double value)
{
  size_t rank1 = findPopIndexRank(id1);
  size_t rank2 = findPopIndexRank(id2);
  size_t rank3 = findPopIndexRank(id3);
  size_t rank4 = findPopIndexRank(id4);
  size_t focalMomIndex = findPi2Index(rank1, rank2, rank3, rank4);

  moments_[focalMomIndex].setValue(value);
}

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
      moments_.push_back(Moment("DD_" + asString(*itI) + "_" + asString(*itJ), 0.));
      moments_.push_back(Moment("H_" + asString(*itI) + "_" + asString(*itJ), 0.));

      for(auto itK = std::begin(popIndices_); itK != std::end(popIndices_); ++itK)
      {
        moments_.push_back(Moment("Dz_" + asString(*itI) + "_" + asString(*itJ) + "_" + asString(*itK), 0.));

        for(auto itL = std::begin(popIndices_); itL != std::end(popIndices_); ++itL)
          moments_.push_back(Moment("pi2_" + asString(*itI) + "_" + asString(*itJ) + "_" + asString(*itK) + "_" + asString(*itL), 0.));
      }
    }
  }

  // NOTE this crucially determines the (lexicographical) order of stats in rows of Y and rows/cols of transition matrices
  std::sort(std::begin(moments_), std::end(moments_), [&](int a, int b) { a.getName() > b.getName() });
}
