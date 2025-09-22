/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 16/09/2025
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

class Model : public bpp::AbstractParameterAliasable, public bpp::FunctionInterface
{

private:
  std::string name_; // model label
  std::vector<std::string> frozenParams_;
  std::vector<std::shared_ptr<Epoch>> epochs_; // each contains its own set of params and operators
  std::shared_ptr<Data> data_;

  MatrixEngine::VectorVariant expected_;
  double compLogLikelihood_;

  // for continuous-time integration
  // (members to make computeExpectedSumStats() work smoothly within fireParameterChanged)
  bool continuousTime_;
  double dt_;
  double totalTime_;
  double errorTolerance_; // error tolerance for adaptive integration

public:
  Model(const std::string& name,
        const std::vector<std::shared_ptr<Epoch>>& epochs,
        std::shared_ptr<Data> data,
        bool continuousTime, double dt, double totTime, double tol):
  AbstractParameterAliasable(""),
  name_(name),
  frozenParams_(),
  epochs_(epochs),
  data_(data),
  expected_(),
  compLogLikelihood_(-1.),
  continuousTime_(continuousTime),
  dt_(dt),
  totalTime_(totTime),
  errorTolerance_(tol)
  {
    for(auto it = std::begin(epochs); it != std::end(epochs); ++it)
      addParameters_((*it)->getParameters());

    linkMoments_();

    if(dt_ <= 0.0 || dt_ > 1.)
      throw bpp::Exception("Model::Invalid time step dt_! Must be in (0, 1].");

    if(totalTime_ <= 0.0 || !std::isfinite(totalTime_))
      throw bpp::Exception("Model::Invalid totalTime_! Must be positive and finite.");

    if(errorTolerance_ <= 0.0 || errorTolerance_ > 1.0)
      throw bpp::Exception("Model::Invalid errorTolerance_! Must be in (0, 1].");
  }

  ~Model()
  {
    std::vector<std::string> paramNames(0);
    paramNames.reserve(getParameters().size());

    for(size_t i = 0; i < getParameters().size(); ++i)
      paramNames.emplace_back(getParameters()[i].getName());

    deleteParameters_(paramNames);
  }

  Model* clone() const override
  {
    return new Model(*this);
  }

  void fireParameterChanged(const bpp::ParameterList& params);

  void setParameters(const bpp::ParameterList& params)
  {
    AbstractParameterAliasable::setParametersValues(params);
  }

  double getValue() const override
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

  std::shared_ptr<Data> getData() const
  {
    return data_;
  }

  bool continuousTime() const
  {
    return continuousTime_;
  }

  double getDt() const
  {
    return dt_;
  }

  double getTotalTime() const
  {
    return totalTime_;
  }

  double getErrorTolerance() const
  {
    return errorTolerance_;
  }

  const MatrixEngine::VectorVariant& getExpectedStats() const
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

  bpp::ParameterList getUnfrozenParameters() const
  {
    bpp::ParameterList unfrozen = getIndependentParameters();

    if(frozenParams_.size() > 0)
    {
      for(auto it = std::begin(frozenParams_); it != std::end(frozenParams_); ++it)
        unfrozen.deleteParameter(*it);
    }

    return unfrozen;
  }

  void computeExpectedSumStats()
  {
    if(continuousTime_)
      computeExpectedSumStatsAdaptive(); // default to adaptive scheme to determine optimal dt

    else
      computeExpectedSumStatsDiscrete();
  }

  void computeExpectedSumStatsDiscrete();

  void computeExpectedSumStatsContinuous();

  void computeExpectedSumStatsAdaptive();

  void printAliasedMomentsPerEpoch(const std::string& modelName) const;

  void printMomentsIntermediate(const std::string& modelName, size_t interval, const std::vector<std::string>& momNames) const;

  void printAliasedMoments(std::ostream& stream);

  void compressParameters(bool aliasOverEpochs, bool aliasOverPops);

private:
  void linkMoments_();

  void updateEpochs_(const bpp::ParameterList& params);

  void computeCompositeLogLikelihood_();
};

#endif
