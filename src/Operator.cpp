/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 05/08/2022
 *
 */


#include <cmath>

#include <Bpp/Numeric/Constraints.h>
#include <Bpp/Numeric/Parameter.h>

#include "Operator.h"


void Operator::fireParameterChanged(const bpp::ParameterList& params)
{
  matchParametersValues(params);
  setUpMatrix();
}
