/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 05/08/2022
 *
 */


#ifndef _MODEL_H_
#define _MODEL_H_

#include <iostream>
#include <cmath>
#include <cstring>
#include <vector>
#include <algorithm>
#include <cstdlib>

#include <Bpp/Numeric/AbstractParameterAliasable.h>
#include <Bpp/Numeric/Constraints.h>
#include <Bpp/Numeric/ParameterList.h>
#include <Bpp/Numeric/Function/Functions.h>

class Model:
  public bpp::AbstractParameterAliasable,
  public bpp::Function
{

private:
  Eigen::SparseMatrix<int, Dynamic, Dynamic> drift_;
  Eigen::SparseMatrix<int, Dynamic, Dynamic> migration_;
  Eigen::SparseMatrix<int, Dynamic, Dynamic> recombination_;
  Eigen::SparseMatrix<int, Dynamic, Dynamic> mutation_;
  Eigen::SparseMatrix<int, Dynamic, Dynamic> selection_;
  Eigen::Matrix<double, Dynamic, Dynamic> combinedOperator_;

  bpp::ParameterList params_;
  std::vector<double> expectedSumStats_;

  double logLikelihood_;
  double aic_;
  
public:
  Model(const bpp::ParameterList& params):
  drift_(),
  migration_(),
  recombination_(),
  mutation_(),
  selection_(),
  combinedOperator_(),
  params_(),
  expectedSumStats_(),
  logLikelihood_(-1.),
  aic_(-1.)
  {
    includeParameters_(params);
  }
  
  Model():
  drift_(),
  migration_(),
  recombination_(),
  mutation_(),
  selection_(),
  combinedOperator_(),
  params_(),
  expectedSumStats_(),
  logLikelihood_(-1.),
  aic_(-1.)
  { }

  Model* clone() const
  {
    return new Model(*this);
  }

  void setParameters(const bpp::ParameterList& params)
  {
    if()
      includeParameters_(params);

    else
      Model::setParametersValues(params);
  }

  double getValue() const
  {
    return -logLikelihood_;
  }
  
  void fireParameterChanged(const bpp::ParameterList& params); // sets updated values
  
  double getLogLikelihood()
  {
    return logLikelihood_;
  }
  
  void computeAic()
  {
    aic_ = 2. * getNumberOfIndependentParameters() - 2. * logLikelihood_;
  }

  double getAic()
  {
    return aic_;
  }

  void update(const bpp::ParameterList& params); // updates matrices based on params

  void computeExpectedSumStats(const AbstractType& data);

  const std::vector<double>& getExpectedSumStats() const
  {
    return expectedSumStats_;
  }

};

#endif




