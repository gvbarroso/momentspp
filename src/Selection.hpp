/*
 * Authors: Gustavo V. Barroso
 * Created: 22/08/2022
 * Last modified: 14/06/2022
 *
 */


#ifndef _SELECTION_H_
#define _SELECTION_H_

#include "AbstractOperator.hpp"
#include "SumStatsLibrary.hpp"

class Selection: public AbstractOperator
{

public:
  Selection(const bpp::ParameterList selParams, const SumStatsLibrary& sslib):
  AbstractOperator(sslib.getPopIndices())
  {
    includeParameters_(selParams);
    prevParams_.addParameters(getParameters()); // inits list of "previous" parameters
    setUpMatrices_(sslib);
  }

  Selection(double val, std::shared_ptr<bpp::IntervalConstraint> ic, const SumStatsLibrary& sslib):
  AbstractOperator(sslib.getPopIndices())
  {
    addParameter_(new bpp::Parameter("s", val, ic));
    prevParams_.addParameters(getParameters()); // inits list of "previous" parameters
    setUpMatrices_(sslib);
  }

  virtual Selection* clone() const override
  {
    return new Selection(*this);
  }

  void setUpMatrices_(const SumStatsLibrary& sslib) override;

  void updateMatrices_() override;

};

#endif
