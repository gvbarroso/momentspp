/*
 * Authors: Gustavo V. Barroso
 * Created: 09/08/2022
 * Last modified: 08/09/2022
 *
 */


#ifndef _RECOMBINATION_H_
#define _RECOMBINATION_H_

#include "Operator.hpp"

class Recombination:
  public Operator
{

public:
  Recombination(const bpp::ParameterList& params, const SumStatsLibrary& ssl, size_t exponent):
  Operator()
  {
    includeParameters_(params);
    prevParams_.addParameters(getParameters()); // inits list of "previous" parameters
    setExponent(exponent);
    setUpMatrices_(ssl);
  }

  void setUpMatrices_(const SumStatsLibrary& ssl);

  void updateMatrices_();

};

#endif
