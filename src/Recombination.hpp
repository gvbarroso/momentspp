/*
 * Authors: Gustavo V. Barroso
 * Created: 09/08/2022
 * Last modified: 10/08/2022
 *
 */


#ifndef _RECOMBINATION_H_
#define _RECOMBINATION_H_

#include "Operator.hpp"

class Recombination:
  public Operator
{

public:
  Recombination():
  Operator()
  { }

  Recombination(const bpp::ParameterList& params):
  Operator(params)
  { }

  Recombination(const bpp::ParameterList& params, size_t matrixSize):
  Operator(params, matrixSize)
  { }

};

#endif
