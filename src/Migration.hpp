/*
 * Authors: Gustavo V. Barroso
 * Created: 10/08/2022
 * Last modified: 22/08/2022
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

  Migration(const bpp::ParameterList& params, const SumStatsLibrary& ssl):
  Operator(params)
  {
    setUpMatrices_(ssl);
  }

};

#endif
