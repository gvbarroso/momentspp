/*
 * Authors: Gustavo V. Barroso
 * Created: 30/08/2022
 * Last modified: 31/08/2022
 *
 */


#ifndef _EPOCH_H_
#define _EPOCH_H_

#include <vector>

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
  // each operator contains bpp parameters and Eigen (sparse) matrices
  std::vector<Operator*> operators_;

  size_t startGen_; // we let the deepest point in relevant time be generation "0"
  size_t endGen_;

public:
  Epoch(const std::vector<Operator*>& operators, size_t start, size_t end, const std::string& name):
  operators_(operators),
  startGen_(sslib),
  endGen_()
  {
    AbstractParameterAliasable::setNamespace(name); // e_0, e_1, e_1 ...
    for(auto it = std::begin(operators); it != std::end(operators); ++it)
      includeParameters_((*it)->getParameters());
  }

  ~Epoch()
  {
    delete ptrs
  }

  Epoch* clone() const
  {
    return new Epoch(*this);
  }

  void fireParameterChanged(const bpp::ParameterList& params);

  void setParameters(const bpp::ParameterList& params)
  {
    AbstractParameterAliasable::setParametersValues(params);

    for(auto it = std::begin(operators_); it != std::end(operators_); ++it)
      (*it)->fireParametersChanged(params);
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

  void computeExpectedSumStats(const Eigen::Matrix<double, Dynamic, Dynamic>& matrix, Eigen::Matrix<double, Dynamic, 1>& y);

private:
  void updateOperators_(const bpp::ParameterList& params);

  Eigen::Matrix<double, Dynamic, Dynamic> integrateOperators_();

};

#endif
