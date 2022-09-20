/*
 * Authors: Gustavo V. Barroso
 * Created: 30/08/2022
 * Last modified: 19/09/2022
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

#include <Bpp/Numeric/AbstractParameterAliasable.h>

#include "AbstractOperator.hpp"
#include "Moment.hpp"
#include "SumStatsLibrary.hpp"

class Epoch:
  public bpp::AbstractParameterAliasable
{

private:
  SumStatsLibrary ssl_;

  // each operator contains Eigen matrices and a subset of the parameters
  std::vector<std::shared_ptr<AbstractOperator>> operators_;
  std::map<size_t, std::shared_ptr<Population>> pops_; // population id -> class object (containing that same id as a member)

  EigenDecomposition eigenDec_;
  Eigen::MatrixXd transitionMatrix_; // all sparse operators combined into a dense matrix
  Eigen::VectorXd steadYstate_; // based on the parameters of this epoch, in practice only the one from deepest epoch will be used as model seed

  size_t startGen_; // we let the deepest point in relevant time be generation "0"
  size_t endGen_;

public:
  Epoch(const SumStatsLibrary& ssl,
        const std::vector<std::shared_ptr<AbstractOperator>>& ops,
        const std::map<size_t, std::shared_ptr<Population>>& pops,
        size_t start, size_t end, const std::string& name):
  bpp::AbstractParameterAliasable(name), // set namespace
  ssl_(ssl),
  operators_(ops),
  pops_(pops),
  eigenDec_(),
  transitionMatrix_(),
  startGen_(start),
  endGen_(end)
  {
    for(auto it = std::begin(ops); it != std::end(ops); ++it)
      includeParameters_((*it)->getParameters()); // NOTE shareParameters

    computeSteadyState_();
  }

  ~Epoch()
  {
    std::cout << "Destruction of Epoch with parameters:\n";
    getParameters().printParameters(std::cout);

    deleteParameters(getParameterNames()); // NOTE does this free memory?
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

  const EigenDecomposition& getEigenDecompositions()
  {
    return eigenDec_;
  }

  void computeExpectedSumStats(Eigen::VectorXd& y)
  {
    transitionMatrix_ * y;
  }

  const SumStatsLibrary& getSslib() const
  {
    return ssl_;
  }

  const std::vector<Moment>& getMoments() const
  {
    return ssl_.getMoments();
  }

private:
  void computeSteadyState_();

  void updateOperators_(const bpp::ParameterList& params)
  {
    for(auto it = std::begin(operators_); it != std::end(operators_); ++it)
      (*it)->fireParameterChanged(params);
  }

};

#endif
