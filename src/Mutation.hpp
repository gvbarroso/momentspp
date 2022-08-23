/*
 * Authors: Gustavo V. Barroso
 * Created: 10/08/2022
 * Last modified: 22/08/2022
 *
 */


#ifndef _MUTATION_H_
#define _MUTATION_H_

#include "Operator.hpp"

class Mutation:
  public Operator
{

public:
  Mutation():
  Operator()
  { }

  Mutation(const bpp::ParameterList& params, const SumStatsLibrary& ssl):
  Operator(params)
  {
    setUpMatrices_(ssl);
  }

};

#endif
