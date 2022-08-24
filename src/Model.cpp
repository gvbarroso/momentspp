/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 24/08/2022
 *
 */


#include "Model.hpp"

void Model::fireParameterChanged(const bpp::ParameterList& params)
{
  matchParametersValues(params);
  updateOperators_(params);
  computeExpectedSumStats_();

  computeCompositeLogLikelihood_(sslib_->getYvec());
  computeAic_();
}

void Model::updateOperators_(const bpp::ParameterList& params)
{
  if(params.getCommonParametersWith(drift_->getIndependentParameters()).size() > 0)
    drift_->fireParameterChanged(params);

  if(params.getCommonParametersWith(migration_->getIndependentParameters()).size() > 0)
    migration_->fireParameterChanged(params);

  if(params.getCommonParametersWith(mutation_->getIndependentParameters()).size() > 0)
    mutation_->fireParameterChanged(params);

  if(params.getCommonParametersWith(recombination_->getIndependentParameters()).size() > 0)
    recombination_->fireParameterChanged(params);

  integrateOperators_(); // updates combinedOperator_
}

void Model::integrateOperators_()
{
  if(continuousTime_) // we combine operators by matrix addition
  {
    combinedOperator_ = drift->getCombinedMatrix() +
                        migration->getCombinedMatrix() +
                        recombination->getCombinedMatrix() +
                        mutation->getCombinedMatrix();
  }

  else // we combine operators by matrix multiplication
  {
    combinedOperator_ = mutation->getCombinedMatrix() *
                        recombination->getCombinedMatrix() *
                        drift->getCombinedMatrix() *
                        migration->getCombinedMatrix();
  }
}

void Model::computeCompositeLogLikelihood_(const Eigen::Matrix<double, Dynamic, 1>& observed)
{
  // observed - expected
}

void Model::computeExpectedSumStats_()
{

}

