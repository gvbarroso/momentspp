/*
 * Authors: Gustavo V. Barroso
 * Created: 10/08/2022
 * Last modified: 10/08/2022
 *
 */


#ifndef _MIGRATION_H_
#define _MIGRATION_H_

#include "Operator.hpp"

class Migration:
  public Operator
{

public:
  Migration():
  Operator()
  { }

  Migration(const bpp::ParameterList& params):
  Operator(params)
  { }

  Migration(const bpp::ParameterList& params, size_t matrixSize):
  Operator(params, matrixSize)
  { }

};

#endif
