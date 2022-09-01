/*
 * Authors: Gustavo V. Barroso
 * Created: 10/08/2022
 * Last modified: 01/09/2022
 *
 */


#ifndef _MUTATION_H_
#define _MUTATION_H_

#include "Operator.hpp"

class Mutation:
  public Operator
{

public:
  Mutation(const bpp::ParameterList& params, const SumStatsLibrary& ssl):
  Operator(params)
  {
    includeParameters_(params);
    prevParams_.addParameters(getParameters()); // inits list of "previous" parameters
    setUpMatrices_(ssl);
  }

};

#endif
