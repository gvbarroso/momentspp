/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 23/08/2022
 *
 */


#include <map>

#include "Model.hpp"

void Model::fireParameterChanged(const bpp::ParameterList& params)
{
  update_(params);
  computeExpectedSumStats_();

  // updates logLikelihood_ and aic_
  computeCompositeLogLikelihood_(expected, observed);
  computeAic_();
}

void Model::update_(const bpp::ParameterList& params)
{

  integrateOperators_();
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

void Model::computeCompositeLogLikelihood_(const std::map<std::string, double>& expected, const std::map<std::string, double>& observed)
{

}

