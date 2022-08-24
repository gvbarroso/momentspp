/*
 * Authors: Gustavo V. Barroso
 * Created: 22/08/2022
 * Last modified: 22/08/2022
 *
 */


#ifndef _SELECTION_H_
#define _SELECTION_H_

#include "Operator.hpp"
#include "SumStatsLibrary.hpp"

class Selection:
  public Operator
{

public:
  Selection():
  Operator()
  { }

  Selection(const bpp::ParameterList& params, const SumStatsLibrary& ssl):
  Operator(params)
  {
    setUpMatrices_(ssl);
  }

};

#endif
