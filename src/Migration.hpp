/*
 * Authors: Gustavo V. Barroso
 * Created: 10/08/2022
 * Last modified: 07/09/2022
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
  Operator()
  {
    // NOTE the constraint that individual migration rates are "small" guaranteed that the rows
    // of the matrix (m_ij's) sum to 1, with main diagonal entries = 1 - sum of values < 1e=5
    includeParameters_(params);
    prevParams_.addParameters(getParameters()); // inits list of "previous" parameters
    setUpMatrices_(ssl);
  }

};

#endif
