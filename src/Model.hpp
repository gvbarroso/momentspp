/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 24/08/2022
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

#include "Drift.hpp"
#include "Migration.hpp"
#include "Mutation.hpp"
#include "Recombination.hpp"
//#include "Selection.hpp"
//#include "Admixture.hpp"
#include "SumStatsLibrary.hpp"

class Model:
  public bpp::AbstractParameterAliasable,
  public bpp::Function
{

private:
  // The fundamental operators derive from base class Operator
  // Each contains matrices of type Eigen::SparseMatrix<double>
  Drift* drift_;
  Migration* migration_;
  Recombination* recombination_;
  Mutation* mutation_;
  //Selection* selection_;
  //Admixture admixture_;

  Eigen::Matrix<double, Dynamic, Dynamic> combinedOperator_; // dense combination of the operators
  Eigen::Matrix<double, Dynamic, 1> expectedY_; // expected sum stats for given parameters

  SumStatsLibrary sslib_;

  // if discrete time, the number of generations spent inside each epoch
  std::vector<size_t> epochGenBoundaries_;

  bool continuousTime_;

  double logLikelihood_;
  double aic_;

public:
  Model(const SumStatsLibrary& sslib):
  drift_(sslib),
  migration_(sslib),
  recombination_(sslib),
  mutation_(sslib),
  //selection_(sslib),
  combinedOperator_(),
  expectedY_(),
  sslib_(sslib),
  epochGenBoundaries_(0),
  continuousTime_(false),
  logLikelihood_(-1.),
  aic_(-1.)
  { }

  Model(const bpp::ParameterList& params, const SumStatsLibrary& sslib):
  drift_(params, sslib),
  migration_(params, sslib),
  recombination_(params, sslib),
  mutation_(params, sslib),
  //selection_(params, sslib),
  combinedOperator_(),
  expectedY_(),
  sslib_(sslib),
  epochGenBoundaries_(0),
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

  void fireParameterChanged(const bpp::ParameterList& params);

  void setParameters(const bpp::ParameterList& params)
  {
    AbstractParameterAliasable::setParametersValues(params);

    drift_->fireParametersChanged(params);
    migration_->fireParametersChanged(params);
    recombination_->fireParametersChanged(params);
    mutation_->fireParametersChanged(params);
    selection_->fireParametersChanged(params);
    admixture_->fireParametersChanged(params);
  }

  double getValue() const
  {
    return -logLikelihood_;
  }

  bool continuousTime()
  {
    return continuousTime_;
  }
  
  double logLikelihood()
  {
    return logLikelihood_;
  }

  double aic()
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

  const Eigen::Matrix<double, Dynamic, 1>& getExpectedSumStats() const
  {
    return expectedY_;
  }

private:
  void integrateOperators_();

  void updateOperators_(const bpp::ParameterList& params);

  void computeExpectedSumStats_(const SomeAbstractType& data);

  void computeCompositeLogLikelihood_(const Eigen::Matrix<double, Dynamic, 1>& observed);

  void computeAic_()
  {
    aic_ = 2. * getNumberOfIndependentParameters() - 2. * logLikelihood_;
  }

};

#endif




