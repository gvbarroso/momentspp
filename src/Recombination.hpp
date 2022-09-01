/*
 * Authors: Gustavo V. Barroso
 * Created: 09/08/2022
 * Last modified: 01/09/2022
 *
 */


#ifndef _RECOMBINATION_H_
#define _RECOMBINATION_H_

#include "Operator.hpp"

class Recombination:
  public Operator
{

public:
  Recombination(const bpp::ParameterList& params, const SumStatsLibrary& ssl):
  Operator()
  {
    includeParameters_(params);
    prevParams_.addParameters(getParameters()); // inits list of "previous" parameters
    setUpMatrices_(ssl);
  }

};

#endif
