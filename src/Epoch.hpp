/*
 * Authors: Gustavo V. Barroso
 * Created: 30/08/2022
 * Last modified: 12/09/2022
 *
 */


#ifndef _EPOCH_H_
#define _EPOCH_H_

#include <vector>
#include <memory>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <Bpp/Numeric/AbstractParameterAliasable.h>

#include "Operator.hpp"

class Epoch:
  public bpp::AbstractParameterAliasable
{

private:
  // each operator contains Eigen matrices and a subset of the parameters
  std::vector<std::shared_ptr<Operator>> operators_;

  EigenDecomposition eigenDec_;
  Eigen::MatrixXd transitionMatrix_;

  size_t startGen_; // we let the deepest point in relevant time be generation "0"
  size_t endGen_;

public:
  Epoch(const std::vector<std::shared_ptr<Operator>>& operators, size_t start, size_t end, const std::string& name):
  bpp::AbstractParameterAliasable(name),
  operators_(operators),
  eigenDec_(),
  transitionMatrix_(),
  startGen_(start),
  endGen_(end)
  {
    for(auto it = std::begin(operators); it != std::end(operators); ++it)
      includeParameters_((*it)->getParameters());
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

private:
  void updateOperators_(const bpp::ParameterList& params)
  {
    for(auto it = std::begin(operators_); it != std::end(operators_); ++it)
      (*it)->fireParameterChanged(params);
  }

};

#endif
