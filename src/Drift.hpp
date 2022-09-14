/*
 * Authors: Gustavo V. Barroso
 * Created: 09/08/2022
 * Last modified: 14/09/2022
 *
 */


#ifndef _DRIFT_H_
#define _DRIFT_H_

#include "Operator.hpp"
#include "SumStatsLibrary.hpp"

class Drift:
  public Operator
{

public:
  Drift(const SumStatsLibrary& ssl):
  Operator()
  {
    // for each population modeled in the epoch this operator belongs to, add Ne parameter
    for(auto itI = std::begin(ssl.getPopMap()); itI != std::end(ssl.getPopMap()); ++itI)
      addParameter(new bpp::Parameter("N_" + bpp::TextTools::toString((*itI).first), 1e+4, bpp::Parameter::R_PLUS_STAR)); // WARNING set this to be >= 1e+3?

    prevParams_.addParameters(getParameters()); // inits list of "previous" parameters
    setUpMatrices_(ssl);
  }

  void setUpMatrices_(const SumStatsLibrary& ssl);

  void updateMatrices_();
};

#endif
