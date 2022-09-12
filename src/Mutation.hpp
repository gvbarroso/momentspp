/*
 * Authors: Gustavo V. Barroso
 * Created: 10/08/2022
 * Last modified: 09/09/2022
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
  Operator()
  {
    includeParameters_(params);
    prevParams_.addParameters(getParameters()); // inits list of "previous" parameters
    setUpMatrices_(ssl);
  }

  void setUpMatrices_(const SumStatsLibrary& ssl);

  void updateMatrices_();

};

#endif
