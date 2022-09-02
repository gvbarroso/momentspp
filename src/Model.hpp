/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 02/09/2022
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

#include "Epoch.h"
#include "SumStatsLibrary.hpp"

class Model:
  public bpp::AbstractParameterAliasable,
  public bpp::Function
{

private:
  std::string name_; // model id
  std::vector<std::shared_ptr<Epoch>> epochs_; // each epoch contains its own parameters and operators
  SumStatsLibrary sslib_; // "Utils" class
  std::vector<std::string> frozenParams_;

  double compLogLikelihood_;

public:
  Model(const std::string& name, const std::vector<std::shared_ptr<Epoch>>& epochs, const SumStatsLibrary& sslib):
  name_(name),
  epochs_(epochs),
  sslib_(sslib),
  frozenParams_(0),
  compLogLikelihood_(-1.)
  {
    for(auto it = std::begin(epochs); it != std::end(epoch); ++it)
      includeParameters_((*it)->getParameters());
  }

  ~Model()
  { }

  Model* clone() const
  {
    return new Model(*this);
  }

  void fireParameterChanged(const bpp::ParameterList& params);

  void setParameters(const bpp::ParameterList& params)
  {
    AbstractParameterAliasable::setParametersValues(params);

    for(auto it = std::begin(epochs_); it != std::end(epochs_); ++it)
      (*it)->fireParametersChanged(params);
  }

  double getValue() const
  {
    return -compLogLikelihood_;
  }
  
  double comLogLikelihood()
  {
    return compLogLikelihood_;
  }

  const std::string& getName()
  {
    return name_;
  }

  const std::vector<Epoch*>& getEpochs()
  {
    return epochs_;
  }

  Eigen::Matrix<double, Dynamic, Dynamic> getCombinedOperator()
  {
    return combinedOperator_;
  }

  void freezeParameter(const std::string& name)
  {
    if(hasParameter(name))
      frozenParams_.push_back(name);

    else
      throw bpp::Exception("Model::Attempted to freeze non-existing parameter " + name);
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
  void updateEpochs_(const bpp::ParameterList& params);

  void computeExpectedSumStats_(const Eigen::Matrix<double, Dynamic, Dynamic>& matrix);

  void computeCompositeLogLikelihood_(const Eigen::Matrix<double, Dynamic, 1>& observed);

};

#endif




