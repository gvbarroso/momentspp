/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 29/07/2022
 *
 */


#include "Model.h"

void Model::fireParameterChanged(const bpp::ParameterList& params)
{
  drift_->update(params);
  migration_->update(params);
  recombination_->update(params);
  mutation_->update(params);
  selection_->update(params);

  combinedOperator_ = combineOperators_();
  combinedOperator_.computeExpectedSumStats(data);
  expectedSumStats_ = combinedOperator_.getExpectedSumStats();

  logLikelihood_ =  fetchCompositeLogLikelihood(expected, observed);
  
  computeAic();
}

