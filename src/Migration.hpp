/*
 * Authors: Gustavo V. Barroso
 * Created: 10/08/2022
 * Last modified:01/09/2022
 *
 */


#ifndef _MIGRATION_H_
#define _MIGRATION_H_

#include "Operator.hpp"

class Migration:
  public Operator
{

public:
  Migration(const bpp::ParameterList& params, const SumStatsLibrary& ssl):
  Operator(params)
  {
    includeParameters_(params);
    prevParams_.addParameters(getParameters()); // inits list of "previous" parameters
    setUpMatrices_(ssl);
  }

};

#endif
