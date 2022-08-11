/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 09/08/2022
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
  // the fundamental operators derive from base class Operator
  // and contain matrices of type Eigen::SparseMatrix<int, Dynamic, Dynamic>
  Drift* drift_;
  Migration* migration_;
  Recombination* recombination_;
  Mutation* mutation_;
  Selection* selection_;

  // this is just a combination of the above operators
  Eigen::Matrix<double, Dynamic, Dynamic> combinedOperator_;

  bpp::ParameterList params_;

  SumStatsLibrary ssl_;

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
  ssl_(),
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
  ssl_(),
  logLikelihood_(-1.),
  aic_(-1.)
  { }

  Model* clone() const
  {
    return new Model(*this);
  }

  void fireParameterChanged(const bpp::ParameterList& params); // sets updated values

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
  
  double getLogLikelihood()
  {
    return logLikelihood_;
  }

  Drift* getDriftOperator()
  {
    return drift_;
  }

  Migration* getMigrationOperator()
  {
    return migration_;
  }

  Eigen::SparseMatrix<int, Dynamic, Dynamic> getRecombinationOperator()
  {
    return recombination_;
  }

  Eigen::SparseMatrix<int, Dynamic, Dynamic> getMutationOperator()
  {
    return mutation_;
  }

  Eigen::SparseMatrix<int, Dynamic, Dynamic> getSelectionOperator()
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

  double getAic()
  {
    return aic_;
  }

private:
  void integrateOperators_();

  void update_(const bpp::ParameterList& params); // updates matrices based on params

  void computeExpectedSumStats_(const AbstractType& data);

};

#endif




