/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 25/08/2022
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
#include <Eigen/Eigenvalues>

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
  // each operator contains parameters and matrices of type Eigen::SparseMatrix<double>
  std::vector<Operator*> operators_;

  SumStatsLibrary sslib_;

  DataType data_;

  // if discrete time, the number of generations spent inside each epoch
  std::vector<size_t> epochGenBoundaries_;
  bool continuousTime_;

  double logLikelihood_;
  double aic_;

public:
  Model(const std::vector<Operator*>& operators, const bpp::ParameterList& params, const SumStatsLibrary& sslib):
  operators_(operators),
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

    for(auto it = std::begin(operators_); it != std::end(operators_); ++it)
      (*it)->fireParametersChanged(params);
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

  const std::vector<Operator*>& getOperators()
  {
    return operators_;
  }

  Eigen::Matrix<double, Dynamic, Dynamic> getCombinedOperator()
  {
    return combinedOperator_;
  }

  const Eigen::Matrix<double, Dynamic, 1>& getExpectedSumStats() const
  {
    return expectedY_;
  }

  void freezeParameter(const std::string& paramName)
  {
    // TODO
  }

private:
  Eigen::Matrix<double, Dynamic, Dynamic> integrateOperators_();

  void updateOperators_(const bpp::ParameterList& params);

  void computeExpectedSumStats_(const Eigen::Matrix<double, Dynamic, Dynamic>& matrix);

  void computeCompositeLogLikelihood_(const Eigen::Matrix<double, Dynamic, 1>& observed);

  void computeAic_()
  {
    aic_ = 2. * getNumberOfIndependentParameters() - 2. * logLikelihood_;
  }

};

#endif




