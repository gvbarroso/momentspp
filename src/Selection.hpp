/*
 * Authors: Gustavo V. Barroso
 * Created: 22/08/2022
 * Last modified: 29/08/2022
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

  Selection(const std::vector<mpfr::mpreal>& vals, std::shared_ptr<bpp::IntervalConstraint> ic, const SumStatsLibrary& sslib):
  AbstractOperator(sslib.getPopIndices())
  {
    // for each population modeled in the epoch *this operator belongs to, add s parameter
    for(size_t i = 0; i < popIndices_.size(); ++i)
      addParameter_(new bpp::Parameter("s_" + bpp::TextTools::toString(popIndices_[i]), vals[i], ic));

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
