/*
 * Authors: Gustavo V. Barroso
 * Created: 30/08/2022
 * Last modified: 30/08/2022
 *
 */


#ifndef _EPOCH_H_
#define _EPOCH_H_

#include <vector>

#include <Bpp/Numeric/AbstractParameterAliasable.h>

#include "Operator.hpp"

class Epoch:
  public bpp::AbstractParameterAliasable
{

private:
  // each operator contains bpp parameters and Eigen (sparse) matrices
  std::vector<Operator*> operators_;

  size_t startGen_;
  size_t endGen_;

public:
  Epoch(const std::vector<Operator*>& operators_, const bpp::ParameterList& params, size_t start, size_t end):
  operators_(operators),
  startGen_(sslib),
  endGen_()
  {
    includeParameters_(params);
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

  size_t length()
  {
    return endGen_ - startGen_;
  }

};

#endif
