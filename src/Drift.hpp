/*
 * Authors: Gustavo V. Barroso
 * Created: 09/08/2022
 * Last modified: 21/09/2022
 *
 */


#ifndef _DRIFT_H_
#define _DRIFT_H_

#include "AbstractOperator.hpp"
#include "SumStatsLibrary.hpp"

class Drift:
  public AbstractOperator
{

public:
  Drift(const SumStatsLibrary& ssl):
  AbstractOperator()
  {
    // for each population modeled in the epoch this operator belongs to, add Ne parameter
    for(auto itI = std::begin(ssl.getPopIndices()); itI != std::end(ssl.getPopIndices()); ++itI)
      addParameter_(new bpp::Parameter("N_" + bpp::TextTools::toString((*itI)), 1e+4, bpp::Parameter::R_PLUS_STAR)); // TODO >= 1e+3?

    prevParams_.addParameters(getParameters()); // inits list of "previous" parameters
    setUpMatrices_(ssl);
  }

  void setUpMatrices_(const SumStatsLibrary& ssl);

  void updateMatrices_();
};

#endif
