/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 29/07/2022
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

#include <Bpp/Numeric/Function/Functions.h>

#include "Operator.h"

class Model:
  public bpp::Function
{

private:
  // NOTE "empty" operator given by the indentity matrix with an empty ParameterList?
  Operator driftOperator_;
  Operator migrationOperator_;
  Operator recombinationOperator_;
  Operator mutationOperator_;
  Operator selectionOperator_;
  Operator combinedOperator_;

  bpp::ParameterList params_;
  std::vector<double> expectedSumStats_; // NOTE here or within each Operator?

  double logLikelihood_;
  double aic_;
  
public:
  Model(const bpp::ParameterList& params):
  Operator(),
  logLikelihood_(-1.),
  aic_(-1.)
  {
    includeParameters_(params);
  }
  
  Model(const bpp::ParameterList& params):
  Operator(),
  logLikelihood_(-1.),
  aic_(-1.)
  {
    includeParameters_(params);
  }
  
  Model():
  Operator(mmsmc, numberOfKnots, splinesType),
  logLikelihood_(-1.),
  aic_(-1.)
  { }

  Model* clone() const
  {
    return new Model(*this);
  }

  void setParameters(const bpp::ParameterList& params)
  {
    Model::setParametersValues(params);
  }

  double getValue() const
  {
    return -logLikelihood_;
  }
  
  void fireParameterChanged(const bpp::ParameterList& params); // sets updated value
  
  double getLogLikelihood()
  {
    return logLikelihood_;
  }
  
  void computeAic()
  {
    aic_ = 2. * getNumberOfIndependentParameters() - 2. * logLikelihood_;
  }

  double getAic()
  {
    return aic_;
  }

};

#endif




