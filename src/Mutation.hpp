/*
 * Authors: Gustavo V. Barroso
 * Created: 10/08/2022
 * Last modified: 10/08/2022
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

  Mutation(const bpp::ParameterList& params):
  Operator(params)
  { }

  Mutation(const bpp::ParameterList& params, size_t matrixSize):
  Operator(params, matrixSize)
  { }

};

#endif
