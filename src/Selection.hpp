/*
 * Authors: Gustavo V. Barroso
 * Created: 22/08/2022
 * Last modified: 05/12/2022
 *
 */


#ifndef _SELECTION_H_
#define _SELECTION_H_

#include "AbstractOperator.hpp"
#include "SumStatsLibrary.hpp"

class Selection:
  public AbstractOperator
{

public:
  Selection(const bpp::ParameterList selParams, const SumStatsLibrary& sslib):
  AbstractOperator(sslib.getNumStats())
  {
    includeParameters_(selParams);
    prevParams_.addParameters(getParameters()); // inits list of "previous" parameters
    setUpMatrices_(sslib);
  }

  Selection(std::shared_ptr<bpp::IntervalConstraint> ic, const SumStatsLibrary& ssl):
  AbstractOperator(sslib.getNumStats())
  {
    // includeParameters_(params);
    prevParams_.addParameters(getParameters()); // inits list of "previous" parameters
    setUpMatrices_(ssl);
  }

};

#endif
