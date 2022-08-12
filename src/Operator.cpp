/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 12/08/2022
 *
 */


#include <cmath>

#include <Bpp/Numeric/Constraints.h>
#include <Bpp/Numeric/Parameter.h>

#include "Operator.hpp"


void Operator::fireParameterChanged(const bpp::ParameterList& params)
{
  matchParametersValues(params);
  update_(); // updates matrix_ and prevParam_
}
