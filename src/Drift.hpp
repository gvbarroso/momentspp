/*
 * Authors: Gustavo V. Barroso
 * Created: 09/08/2022
 * Last modified: 09/08/2022
 *
 */


#ifndef _DRIFT_H_
#define _DRIFT_H_

#include "Operator.hpp"

class Drift:
  public Operator
{

public:
  Drift():
  Operator()
  { }

  Drift(const bpp::ParameterList& params):
  Operator(params)
  { }

  Drift(const bpp::ParameterList& params, size_t matrixSize):
  Operator(params, matrixSize)
  { }

};

#endif
