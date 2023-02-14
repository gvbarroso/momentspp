/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 13/02/2023
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
#include <utility>
#include <map>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <Bpp/Numeric/AbstractParameterAliasable.h>
#include <Bpp/Numeric/Constraints.h>
#include <Bpp/Numeric/ParameterList.h>
#include <Bpp/Numeric/Function/Functions.h>

#include "Epoch.hpp"
#include "Data.hpp"

class Model: public bpp::AbstractParameterAliasable, public bpp::Function
{

private:
  std::string name_; // model label / id
  std::vector<std::shared_ptr<Epoch>> epochs_; // each epoch contains its own parameters and operators
  std::vector<std::string> frozenParams_;

  Data data_;
  Eigen::VectorXd expected_;

  double compLogLikelihood_;

public:
  Model(const std::string& name, const std::vector<std::shared_ptr<Epoch>>& epochs, const Data& data):
  AbstractParameterAliasable(""),
  name_(name),
  epochs_(epochs),
  frozenParams_(0),
  data_(data),
  expected_(),
  compLogLikelihood_(-1.)
  {
    for(auto it = std::begin(epochs); it != std::end(epochs); ++it)
      shareParameters_((*it)->getParameters());

    linkMoments_();
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
  }

  double getValue() const
  {
    return -compLogLikelihood_;
  }
  
  double compLogLikelihood()
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

  const Eigen::VectorXd& getExpectedStats()
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

  void unfreezeParameter(const std::string& name)
  {
    if(hasParameter(name))
    {
      auto it = std::find(std::begin(frozenParams_), std::end(frozenParams_), name);

      if(it != std::end(frozenParams_))
        it = frozenParams_.erase(it);

      else
        throw bpp::Exception("Model::Attempted to unfreeze unfrozen parameter " + name);
    }

    else
      throw bpp::Exception("Model::Attempted to unfreeze non-existing parameter " + name);
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

  void computeExpectedSumStats();

  void aliasMoments();

  void printAliasedMoments();

private:
  void linkMoments_();

  void admixture_();

  void updateEpochs_(const bpp::ParameterList& params);

  void computeCompositeLogLikelihood_(const Eigen::VectorXd& obsMeans, const Eigen::MatrixXd& obsCovarMat);

};

#endif




