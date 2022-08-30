/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 30/08/2022
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
  std::string name_; // model id

  // each operator contains bpp parameters and Eigen (sparse) matrices
  std::vector<Operator*> operators_; // NOTE: make this a vector of Epochs?!
  SumStatsLibrary sslib_;
  DataType data_;

  std::vector<std::string> frozenParams_;
  // if discrete time, the number of generations spent inside each epoch
  std::vector<size_t> epochGenBoundaries_;
  bool continuousTime_;

  double compLogLikelihood_;
  double aic_;

public:
  Model(const std::string& name, const std::vector<Operator*>& operators, const bpp::ParameterList& params, const SumStatsLibrary& sslib):
  name_(name)
  operators_(operators),
  sslib_(sslib),
  data_(),
  frozenParams_(0),
  epochGenBoundaries_(0),
  continuousTime_(false),
  compLogLikelihood_(-1.),
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
    return -compLogLikelihood_;
  }

  bool continuousTime()
  {
    return continuousTime_;
  }
  
  double comLogLikelihood()
  {
    return compLogLikelihood_;
  }

  double aic()
  {
    return aic_;
  }

  const std::string& getName()
  {
    return name_;
  }

  const std::vector<Operator*>& getOperators()
  {
    return operators_;
  }

  Eigen::Matrix<double, Dynamic, Dynamic> getCombinedOperator()
  {
    return combinedOperator_;
  }

  void freezeParameter(const std::string& paramName)
  {
    if(hasParameter(paramName))
      frozenParams_.push_back(paramName);

    else
      throw bpp::Exception("Model::Attempted to freeze non-existing parameter " + paramName);
  }

  const bpp::ParameterList& getUnfrozenParameters()
  {
    if(frozenParameters.size() == 0)
      return getIndependentParameters();

    else
    {
      bpp::ParameterList unfrozen = getIndependentParameters();
      for(auto it = std::begin(frozenParameters_); it != std::end(frozenParameters); ++it)
        unfrozen.deleteParameter(*it);

      return unfrozen;
    }
  }

private:
  Eigen::Matrix<double, Dynamic, Dynamic> integrateOperators_();

  void updateOperators_(const bpp::ParameterList& params);

  void computeExpectedSumStats_(const Eigen::Matrix<double, Dynamic, Dynamic>& matrix);

  void computeCompositeLogLikelihood_(const Eigen::Matrix<double, Dynamic, 1>& observed);

  void computeAic_()
  {
    aic_ = 2. * getNumberOfIndependentParameters() - 2. * compLogLikelihood_;
  }

};

#endif




