/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 29/07/2022
 *
 */


#include <cmath>

#include <Bpp/Numeric/Constraints.h>
#include <Bpp/Numeric/Parameter.h>

#include "Operator.h"


void Operator::fireParameterChanged(const bpp::ParameterList& params)
{
  update(params);
}

void Operator::update(const bpp::ParameterList& params)
{
  matchParametersValues(params);
  recomputeEntries_(); // of own matrix
}
