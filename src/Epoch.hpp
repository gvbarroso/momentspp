/*
 * Authors: Gustavo V. Barroso
 * Created: 30/08/2022
 * Last modified: 29/06/2023
 *
 */


#ifndef _EPOCH_H_
#define _EPOCH_H_

#include <iostream>
#include <vector>
#include <memory>
#include <utility>
#include <algorithm>
#include <map>

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <eigen3/unsupported/Eigen/MatrixFunctions>

#include <Bpp/Numeric/AbstractParameterAliasable.h>

#include "AbstractOperator.hpp"
#include "Admixture.hpp"
#include "Mutation.hpp"
#include "SumStatsLibrary.hpp"
#include "Population.hpp"

class Epoch: public bpp::AbstractParameterAliasable
{

private:
  std::string name_;
  SumStatsLibrary ssl_; // *this epoch has its own set of moments using its population indices

  // generations ago, from past to present
  size_t startGen_;
  size_t endGen_;

  std::vector<std::shared_ptr<Population>> pops_;
  std::vector<std::shared_ptr<AbstractOperator>> operators_; // each operator contains matrices and a subset of the parameters

  Eigen::MatrixXd transitionMatrix_; // all sparse operators combined into a dense matrix
  Eigen::VectorXd steadYstate_; // based on the parameters of *this epoch

public:
  Epoch():
  bpp::AbstractParameterAliasable(""),
  name_(),
  ssl_(),
  startGen_(0),
  endGen_(0),
  pops_(0),
  operators_(0),
  transitionMatrix_(),
  steadYstate_()
  { }

  Epoch(const std::string& name, SumStatsLibrary& ssl, size_t start, size_t end,
        const std::vector<std::shared_ptr<AbstractOperator>>& ops,
        const std::vector<std::shared_ptr<Population>>& pops):
  bpp::AbstractParameterAliasable(""),
  name_(name),
  ssl_(ssl),
  startGen_(start),
  endGen_(end),
  pops_(pops),
  operators_(ops),
  transitionMatrix_(),
  steadYstate_()
  {
    for(auto it = std::begin(operators_); it != std::end(operators_); ++it)
      addParameters_((*it)->getParameters());

    bpp::AbstractParameterAliasable::setNamespace(name + ".");
    init_();
  }

  ~Epoch()
  {
    std::vector<std::string> paramNames(0);
    paramNames.reserve(getParameters().size());

    for(size_t i = 0; i < getParameters().size(); ++i)
      paramNames.emplace_back(getParameters()[i].getName());

    deleteParameters_(paramNames);
  }

  Epoch* clone() const
  {
    return new Epoch(*this);
  }

  void fireParameterChanged(const bpp::ParameterList& params);

  void setParameters(const bpp::ParameterList& params)
  {
    bpp::AbstractParameterAliasable::setParametersValues(params);
  }

  const std::string& getName()
  {
    return name_;
  }

  size_t start()
  {
    return startGen_;
  }

  size_t end()
  {
    return endGen_;
  }

  size_t duration()
  {
    return startGen_ - endGen_;
  }

  const Eigen::VectorXd& getSteadyState()
  {
    return steadYstate_;
  }

  const Eigen::MatrixXd& getTransitionMatrix()
  {
    return transitionMatrix_;
  }

  size_t getNumPops()
  {
    return pops_.size();
  }

  const std::vector<std::shared_ptr<Population>>& getPops()
  {
    return pops_;
  }

  void printAttributes(std::ostream& stream)
  {
    stream << name_ << ", from " << startGen_ << " to " << endGen_ << "\n";

    for(auto it = std::begin(pops_); it != std::end(pops_); ++it)
    {
      stream << "\t";
      (*it)->printAttributes(stream);
    }
  }

  std::shared_ptr<Population> fetchPop(size_t id)
  {
    std::shared_ptr<Population> pop = nullptr;
    for(auto it = std::begin(pops_); it != std::end(pops_); ++it)
    {
      if((*it)->getId() == id)
        pop = (*it);
    }

    assert(pop != nullptr);
    return pop;
  }

  std::shared_ptr<Population> fetchPop(const std::string& name)
  {
    std::shared_ptr<Population> pop = nullptr;
    for(auto it = std::begin(pops_); it != std::end(pops_); ++it)
    {
      if((*it)->getName() == name)
        pop = (*it);
    }

    assert(pop != nullptr);
    return pop;
  }

  const SumStatsLibrary& getSslib() const
  {
    return ssl_;
  }

  SumStatsLibrary& getSslib()
  {
    return ssl_;
  }

  const std::vector<std::shared_ptr<Moment>>& getMoments() const
  {
    return ssl_.getMoments();
  }

  std::vector<std::shared_ptr<Moment>>& getMoments()
  {
    return ssl_.getMoments();
  }

  const std::vector<std::shared_ptr<Moment>>& getBasis() const
  {
    return ssl_.getBasis();
  }

  std::vector<std::shared_ptr<Moment>>& getBasis()
  {
    return ssl_.getBasis();
  }

  std::vector<size_t> fetchSelectedPopIds(); // for *this epoch

  void computeExpectedSumStats(Eigen::VectorXd& y);

  void transferStatistics(Eigen::VectorXd& y);

  void updateMoments(const Eigen::VectorXd& y);

  void printRecursions(std::ostream& stream);

  void printTransitionMat(const std::string& fileName) const;

  void pseudoSteadyState();

  void computeSteadyState();

  void testSteadyState();

private:
  void init_();

  void updateOperators_(const bpp::ParameterList& params)
  {
    for(auto it = std::begin(operators_); it != std::end(operators_); ++it)
      (*it)->fireParameterChanged(params);
  }

};

#endif
