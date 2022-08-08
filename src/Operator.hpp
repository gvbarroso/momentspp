/*
 * Authors: Gustavo V. Barroso
 * Created:29/07/2022
 * Last modified: 29/07/2022
 *
 */


#ifndef _OPERATOR_H_
#define _OPERATOR_H_

#include <Bpp/Numeric/Function/Functions.h>
#include <Bpp/Numeric/AbstractParameterAliasable.h>
#include <Bpp/Numeric/Constraints.h>
#include <Bpp/Numeric/ParameterList.h>

class Operator:
  public bpp::AbstractParameterAliasable
{

private:
  Eigen::SparseMatrix<int, Dynamic, Dynamic> matrix_;

public:
  Operator():
  AbstractParameterAliasable("")
  { }

  Operator(const bpp::ParameterList& params):
  AbstractParameterAliasable("")
  {
    bpp::addParameters_(params);

    setUpMatrix();
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

  void update(const bpp::ParameterList& params); // updates matrix based on params

  virtual void setUpMatrix();

};

#endif
