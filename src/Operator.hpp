/*
 * Authors: Gustavo V. Barroso
 * Created:29/07/2022
 * Last modified: 11/08/2022
 *
 */


#ifndef _OPERATOR_H_
#define _OPERATOR_H_


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
  Eigen::SparseMatrix<int, Dynamic, Dynamic> matrix_;

public:
  Operator():
  AbstractParameterAliasable(""),
  matrix_()
  { }

  Operator(const bpp::ParameterList& params):
  AbstractParameterAliasable(""),
  matrix_()
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

  const Eigen::SparseMatrix<int, Dynamic, Dynamic>& getMatrix()
  {
    return matrix_;
  }

  virtual void setUpMatrix(const SumStatsLibrary& sslib);

};

#endif
