/*
 * Authors: Gustavo V. Barroso
 * Created: 09/08/2022
 * Last modified: 22/08/2022
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

  Recombination(const bpp::ParameterList& params, const SumStatsLibrary& ssl):
  Operator(params)
  {
    setUpMatrices_(ssl);
  }

};

#endif
