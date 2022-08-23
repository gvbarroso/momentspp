/*
 * Authors: Gustavo V. Barroso
 * Created: 09/08/2022
 * Last modified: 23/08/2022
 *
 */


#ifndef _DRIFT_H_
#define _DRIFT_H_

#include

#include "Operator.hpp"
#include "SumStatsLibrary.hpp"

class Drift:
  public Operator
{

public:
  Drift(const SumStatsLibrary& ssl):
  Operator()
  {
    for(size_t i = 0; i < ssl->getNumPops(); ++i)
      bpp::addParameter_(new bpp::Parameter("N_" + bpp::toString(i), 1., Parameter::R_PLUS_STAR));

    prevParams_.addParameters(getParameters()); // inits list of "previous" parameters
    setUpMatrices_(ssl);
  }

  Drift(const bpp::ParameterList& params, const SumStatsLibrary& ssl):
  Operator()
  {
    for(size_t i = 0; i < ssl->getNumPops(); ++i)
      bpp::addParameter_(new bpp::Parameter("N_" + bpp::toString(i), 1., Parameter::R_PLUS_STAR));

    matchParametersValues(params);
    prevParams_.addParameters(getParameters()); // inits list of "previous" parameters
    setUpMatrices_(ssl);
  }

};

#endif
