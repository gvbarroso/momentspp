/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 15/09/2025
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

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <eigen3/unsupported/Eigen/MPRealSupport>

#include <Bpp/Numeric/AbstractParameterAliasable.h>
#include <Bpp/Numeric/Constraints.h>
#include <Bpp/Numeric/ParameterList.h>
#include <Bpp/Numeric/Function/Functions.h>

#include "Vector.hpp"
#include "Epoch.hpp"
#include "Data.hpp"

class Model: public bpp::AbstractParameterAliasable, public bpp::FunctionInterface
{

private:
  std::string name_; // model label
  std::vector<std::string> frozenParams_;
  std::vector<std::shared_ptr<Epoch>> epochs_; // each contains its own set of params and operators
  std::shared_ptr<Data> data_;

  MatrixEngine::VectorVariantEigen expected_;
  double compLogLikelihood_;

public:
  Model(const std::string& name, const std::vector<std::shared_ptr<Epoch>>& epochs, std::shared_ptr<Data> data = nullptr):
  AbstractParameterAliasable(""),
  name_(name),
  frozenParams_(0),
  epochs_(epochs),
  data_(data),
  expected_(0),
  compLogLikelihood_(-1.)
  {
    for(auto it = std::begin(epochs); it != std::end(epochs); ++it)
      addParameters_((*it)->getParameters());

    linkMoments_();
  }

  ~Model()
  {
    std::vector<std::string> paramNames(0);
    paramNames.reserve(getParameters().size());

    for(size_t i = 0; i < getParameters().size(); ++i)
      paramNames.emplace_back(getParameters()[i].getName());

    deleteParameters_(paramNames);
  }

  Model* clone() const
  {
    return new Model(*this);
  }

  void fireParameterChanged(const bpp::ParameterList& params);

  void setParameters(const bpp::ParameterList& params)
  {
    AbstractParameterAliasable::setParametersValues(params);
  }

  double getValue() const override // NOTE override added on 15/09/2025
  {
    return -compLogLikelihood_;
  }
  
  const std::string& getName() const
  {
    return name_;
  }

  const std::vector<std::shared_ptr<Epoch>>& getEpochs() const
  {
    return epochs_;
  }

  const MatrixEngine::VectorVariantEigen& getExpectedStats() const
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

  void computeExpectedSumStats(bool continuousTime)
  {
    if(continuousTime)
      computeExpectedSumStatsAdaptive(); // default is to use adaptive scheme to determine optimal dt

    else
      computeExpectedSumStatsDiscrete();
  }

  void computeExpectedSumStatsDiscrete();

  void computeExpectedSumStatsContinuous();

  void computeExpectedSumStatsAdaptive();

  void printAliasedMomentsPerEpoch(const std::string& modelName);

  void printMomentsIntermediate(const std::string& modelName, size_t interval, const std::vector<std::string>& momNames);

  void printAliasedMoments(std::ostream& stream);

  void compressParameters(bool aliasOverEpochs, bool aliasOverPops);

private:
  void linkMoments_();

  void updateEpochs_(const bpp::ParameterList& params);

  void computeCompositeLogLikelihood_();

};

#endif




