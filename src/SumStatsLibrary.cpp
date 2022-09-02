/*
 * Authors: Gustavo V. Barroso
 * Created: 05/08/2022
 * Last modified: 02/09/2022
 *
 */

#include "SumStatsLibrary.hpp"


void SumStatsLibrary::init(const PolymorphismData& dataset)
{
  if(initialized_)
    throw bpp::Exception("SumStatsLibrary::attempted to re-initialized library!");

  else
  {
    setNumPops(dataset.getNumPops());
    setOrder(dataset.getOrder());

    includeHetStats_();
    includeLdStats_();

    y_.resize(stats_.size());
    size_t index = 0;

    // inits y_ from stats_ so that they're in same order (we also init matrices from stats_, see Operator::setUpMatrices_)
    for(auto it = std::begin(stats_); it != std::end(stats_); ++it)
    {
      y_(index, 0) = (*it)->second;
      ++index;
    }

    initialized_ = true;
  }
}

void SumStatsLibrary::includeHetStats_()
{
  for(size_t i = 0; i < numPops_; ++i)
    for(size_t j = 0; j < numPops_; ++j)
      stats_.try_emplace("H_" + asString(i) + asString(j), 0.);
}

void SumStatsLibrary::includeLdStats_()
{
  for(size_t i = 0; i < numPops_; ++i)
    for(size_t j = 0; j < numPops_; ++j)
      stats_.try_emplace("DD_" + asString(i) + asString(j), 0.);

  for(size_t i = 0; i < numPops_; ++i)
    for(size_t j = 0; j < numPops_; ++j)
      for(size_t k = 0; k < numPops_; ++k)
        stats_.try_emplace("Dz_" + asString(i) + asString(j) + asString(k), 0.);

  for(size_t i = 0; i < numPops_; ++i)
    for(size_t j = 0; j < numPops_; ++j)
      for(size_t k = 0; k < numPops_; ++k)
        for(size_t l = 0; l < numPops_; ++l)
          stats_.try_emplace("pi2_" + asString(i) + asString(j) + ";" + asString(k) + asString(l));
}
