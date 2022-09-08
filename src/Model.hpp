/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 08/09/2022
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
#include <memory>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <Bpp/Numeric/AbstractParameterAliasable.h>
#include <Bpp/Numeric/Constraints.h>
#include <Bpp/Numeric/ParameterList.h>
#include <Bpp/Numeric/Function/Functions.h>

#include "Epoch.hpp"
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

  Eigen::Matrix<double, Eigen::Dynamic, 1> steadYstate_;
  Eigen::Matrix<double, Eigen::Dynamic, 1> expected_;

  double compLogLikelihood_;

public:
  Model(const std::string& name, const std::vector<std::shared_ptr<Epoch>>& epochs, const SumStatsLibrary& sslib):
  AbstractParameterAliasable(""),
  name_(name),
  epochs_(epochs),
  sslib_(sslib),
  frozenParams_(0),
  steadYstate_(),
  expected_(),
  compLogLikelihood_(-1.)
  {
    for(auto it = std::begin(epochs); it != std::end(epochs); ++it)
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
      (*it)->fireParameterChanged(params);
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

  const std::vector<std::shared_ptr<Epoch>>& getEpochs()
  {
    return epochs_;
  }

  const Eigen::Matrix<double, Eigen::Dynamic, 1>& getSteadyState()
  {
    return steadYstate_;
  }

  const Eigen::Matrix<double, Eigen::Dynamic, 1>& getExpectedStats()
  {
    return expected_;
  }

  void freezeParameter(const std::string& name)
  {
    if(hasParameter(name))
      frozenParams_.push_back(name);

    else
      throw bpp::Exception("Model::Attempted to freeze non-existing parameter " + name);
  }

  bpp::ParameterList getUnfrozenParameters()
  {
    bpp::ParameterList unfrozen = getIndependentParameters();
    if(frozenParams_.size() > 0)
    {
      for(auto it = std::begin(frozenParams_); it != std::end(frozenParams_); ++it)
        unfrozen.deleteParameter(*it);
    }

    return unfrozen;
  }

  void computeSteadyState();

private:
  void updateEpochs_(const bpp::ParameterList& params);

  void computeExpectedSumStats_();

  void computeCompositeLogLikelihood_(const Eigen::VectorXd& obsMeans, const Eigen::MatrixXd& obsCovarMat);

};

#endif




