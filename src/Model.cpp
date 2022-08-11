/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 05/08/2022
 *
 */


#include "Model.hpp"

void Model::fireParameterChanged(const bpp::ParameterList& params)
{
  update_(params);
  integrateOperators_();
  computeExpectedSumStats_();

  // updates logLikelihood_ and aic_
  computeCompositeLogLikelihood_(expected, observed);
  computeAic_();
}

