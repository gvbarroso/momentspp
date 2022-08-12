/*
 * Authors: Gustavo V. Barroso
 * Created:29/07/2022
 * Last modified: 12/08/2022
 *
 */


#ifndef _OPERATOR_H_
#define _OPERATOR_H_

#include <vector>
#include <string>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Dense>

#include <Bpp/Numeric/Function/Functions.h>
#include <Bpp/Numeric/AbstractParameterAliasable.h>
#include <Bpp/Numeric/Constraints.h>
#include <Bpp/Numeric/ParameterList.h>

#include "SumStatsLibrary.hpp"

class Operator:
  public bpp::AbstractParameterAliasable
{

private:
  std::vector<Eigen::SparseMatrix<double, Dynamic, Dynamic>> matrices_; // one matrix per population
  Eigen::SparseMatrix<double, Dynamic, Dynamic> combinedPopMatrix_;

  std::vector<double> prevParams_; // parameters values in immediately previous iteration of optimization

public:
  Operator():
  AbstractParameterAliasable(""),
  matrix_(0)
  { }

  Operator(const bpp::ParameterList& params):
  AbstractParameterAliasable(""),
  matrix_(0)
  {
    bpp::addParameters_(params);
  }

public:
  Operator* clone() const
  {
    return new Operator(*this);
  }

  void setParameters(const bpp::ParameterList& params)
  {
    AbstractParameterAliasable::setParametersValues(params);
  }

  void fireParameterChanged(const bpp::ParameterList& params);

  const Eigen::SparseMatrix<int, Dynamic, Dynamic>& getPopMatrices()
  {
    return matrices_;
  }

  const Eigen::SparseMatrix<int, Dynamic, Dynamic>& getPopMatrix(size_t popIndex)
  {
    return matrices_[popIndex];
  }

  const Eigen::SparseMatrix<int, Dynamic, Dynamic>& getCombinedPopMatrix()
  {
    return combinedPopMatrix_;
  }

  virtual void setUpMatrices(const SumStatsLibrary& sslib);

private:
  virtual void update_();

};

#endif
