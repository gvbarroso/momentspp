/*
 * Authors: Gustavo V. Barroso
 * Created: 09/08/2022
 * Last modified: 22/08/2022
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
  Drift():
  Operator()
  { }

  Drift(const bpp::ParameterList& params, const SumStatsLibrary& ssl):
  Operator(params)
  {
    setUpMatrices_(ssl);
  }

};

#endif
