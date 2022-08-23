/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 23/08/2022
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

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Dense>

#include <Bpp/Numeric/AbstractParameterAliasable.h>
#include <Bpp/Numeric/Constraints.h>
#include <Bpp/Numeric/ParameterList.h>
#include <Bpp/Numeric/Function/Functions.h>

class Model:
  public bpp::AbstractParameterAliasable,
  public bpp::Function
{

private:
  // The fundamental operators derive from base class Operator
  // Each contains matrices of type Eigen::SparseMatrix<int, Dynamic, Dynamic>
  Drift* drift_;
  Migration* migration_;
  Recombination* recombination_;
  Mutation* mutation_;
  Selection* selection_;

  // this is a handy combination of the above operators
  Eigen::Matrix<double, Dynamic, Dynamic> combinedOperator_;

  SumStatsLibrary sslib_; // contains summary statistics

  bool continuousTime_;

  double logLikelihood_;
  double aic_;

public:
  Model():
  drift_(),
  migration_(),
  recombination_(),
  mutation_(),
  selection_(),
  combinedOperator_(),
  sslib_(),
  continuousTime_(false),
  logLikelihood_(-1.),
  aic_(-1.)
  { }

  Model(const bpp::ParameterList& params, const SumStatsLibrary& sslib):
  drift_(),
  migration_(),
  recombination_(),
  mutation_(),
  selection_(),
  combinedOperator_(),
  sslib_(sslib),
  continuousTime_(false),
  logLikelihood_(-1.),
  aic_(-1.)
  {
    includeParameters_(params);
  }

  Model* clone() const
  {
    return new Model(*this);
  }

  void fireParameterChanged(const bpp::ParameterList& params); // sets updated values

  void setParameters(const bpp::ParameterList& params)
  {
    Model::setParametersValues(params);
  }

  double getValue() const
  {
    return -logLikelihood_;
  }

  bool continuousTime()
  {
    return continuousTime_;
  }
  
  double getLogLikelihood()
  {
    return logLikelihood_;
  }

  double getAic()
  {
    return aic_;
  }

  Drift* getDriftOperator()
  {
    return drift_;
  }

  Migration* getMigrationOperator()
  {
    return migration_;
  }

  Recombination* getRecombinationOperator()
  {
    return recombination_;
  }

  Mutation* getMutationOperator()
  {
    return mutation_;
  }

  Selection* getSelectionOperator()
  {
    return selection_;
  }

  Eigen::Matrix<double, Dynamic, Dynamic> getCombinedOperator()
  {
    return combinedOperator_;
  }

  const std::vector<double>& getExpectedSumStats() const
  {
    return expectedSumStats_;
  }
  
  void computeAic()
  {
    aic_ = 2. * getNumberOfIndependentParameters() - 2. * logLikelihood_;
  }

private:
  void integrateOperators_();

  void update_(const bpp::ParameterList& params); // updates matrices based on params

  void computeExpectedSumStats_(const SomeAbstractType& data);

  void computeCompositeLogLikelihood_(const SomeAbstractType& expected, const SomeAbstractType& observed);

};

#endif




