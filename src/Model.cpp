/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 29/07/2022
 *
 */


#include "Model.h"

void Model::fireParameterChanged(const bpp::ParameterList& params)
{
  driftOperator_->update(params);
  migrationOperator_->update(params);
  recombinationOperator_->update(params);
  mutationOperator_->update(params);
  selectionOperator_->update(params);

  combinedOperator_ = combineOperators_();
  combinedOperator_.computeExpectedSumStats(data);
  expectedSumStats_ = combinedOperator_

  logLikelihood_ =  fetchCompositeLogLikelihood(expected, observed);
  
  computeAic();
}

