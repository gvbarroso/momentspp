/*
 * Authors: Gustavo V. Barroso
 * Created: 30/08/2022
 * Last modified: 05/12/2022
 *
 */


#ifndef _EPOCH_H_
#define _EPOCH_H_

#include <vector>
#include <memory>
#include <utility>
#include <algorithm>
#include <map>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/MatrixFunctions> // for es_.pseudoEigenvalueMatrix().pow(duration())

#include <Bpp/Numeric/AbstractParameterAliasable.h>

#include "AbstractOperator.hpp"
#include "Mutation.hpp"
#include "SumStatsLibrary.hpp"
#include "Population.hpp"

class Epoch:
  public bpp::AbstractParameterAliasable
{

private:
  SumStatsLibrary ssl_; // *this epoch has its own set of moments following *this population indices

  size_t startGen_; // we let the deepest point in relevant time be generation '0'
  size_t endGen_;

  // each operator contains Eigen matrices and a subset of the parameters
  std::vector<std::shared_ptr<AbstractOperator>> operators_;
  std::map<size_t, std::shared_ptr<Population>> pops_; // pop-id->class object (containing that same id as a member variable)

  Eigen::MatrixXd transitionMatrix_; // all sparse operators combined into a dense matrix
  Eigen::VectorXd steadYstate_; // based on the parameters of *this epoch

public:
  Epoch():
  bpp::AbstractParameterAliasable(""),
  ssl_(),
  startGen_(0),
  endGen_(0),
  operators_(0),
  pops_(),
  transitionMatrix_(),
  steadYstate_()
  { }

  Epoch(const SumStatsLibrary& ssl, size_t start, size_t end, const std::string& name,
        const std::vector<std::shared_ptr<AbstractOperator>>& ops,
        const std::map<size_t, std::shared_ptr<Population>>& pops):
  bpp::AbstractParameterAliasable(""), // set namespace TODO use name somehow
  ssl_(ssl),
  startGen_(start),
  endGen_(end),
  operators_(ops),
  pops_(pops),
  transitionMatrix_(),
  steadYstate_()
  {
    for(auto it = std::begin(operators_); it != std::end(operators_); ++it)
      shareParameters_((*it)->getParameters());

    computeSteadyState_(); // updates moments inside ssl_
  }

  ~Epoch()
  { }

  Epoch* clone() const
  {
    return new Epoch(*this);
  }

  void fireParameterChanged(const bpp::ParameterList& params);

  void setParameters(const bpp::ParameterList& params)
  {
    bpp::AbstractParameterAliasable::setParametersValues(params);
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
    return endGen_ - startGen_;
  }

  const Eigen::VectorXd& getSteadyState()
  {
    return steadYstate_;
  }

  const Eigen::MatrixXd& getTransitionMatrix()
  {
    return transitionMatrix_;
  }

  const std::map<size_t, std::shared_ptr<Population>>& getPops()
  {
    return pops_;
  }

  const SumStatsLibrary& getSslib() const
  {
    return ssl_;
  }

  SumStatsLibrary& getSslib()
  {
    return ssl_;
  }

  const std::vector<Moment>& getMoments() const
  {
    return ssl_.getMoments();
  }

  std::vector<Moment>& getMoments()
  {
    return ssl_.getMoments();
  }

  void computeExpectedSumStats(Eigen::VectorXd& y);

  void transferStatistics(Eigen::VectorXd& y);

  void updateMoments(const Eigen::VectorXd& y);

  void timeTest(size_t g);

private:
  void computeSteadyState_();

  void updateOperators_(const bpp::ParameterList& params)
  {
    for(auto it = std::begin(operators_); it != std::end(operators_); ++it)
      (*it)->fireParameterChanged(params);
  }

};

#endif
