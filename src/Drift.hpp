/*
 * Authors: Gustavo V. Barroso
 * Created: 09/08/2022
 * Last modified: 08/09/2022
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
  Drift(const bpp::ParameterList& params, const SumStatsLibrary& ssl, size_t exponent):
  Operator()
  {
    includeParameters_(params);
    prevParams_.addParameters(getParameters()); // inits list of "previous" parameters
    setUpMatrices_(ssl);

    setExponent(exponent);
  }

  void setUpMatrices_(const SumStatsLibrary& ssl);

  void updateMatrices_();
};

#endif
