/*
 * Authors: Gustavo V. Barroso
 * Created: 09/08/2022
 * Last modified: 20/04/2023
 *
 */


#ifndef _DRIFT_H_
#define _DRIFT_H_

#include "AbstractOperator.hpp"
#include "SumStatsLibrary.hpp"

class Drift: public AbstractOperator
{

public:
  Drift(const bpp::ParameterList driftParams, const SumStatsLibrary& sslib):
  AbstractOperator()
  {
    includeParameters_(driftParams);
    prevParams_.addParameters(getParameters()); // inits list of "previous" parameters
    setUpMatrices_(sslib);
  }

  Drift(const std::vector<double>& vals, std::shared_ptr<bpp::IntervalConstraint> ic, const SumStatsLibrary& sslib):
  AbstractOperator()
  {
    // for each population modeled in the epoch *this operator belongs to, add Ne parameter
    for(size_t i = 0; i < sslib.getPopIndices().size(); ++i)  // i-th coal rate corresponds top popIndex[i]
      addParameter_(new bpp::Parameter("1/2N_" + bpp::TextTools::toString(i), vals[i], ic));

    prevParams_.addParameters(getParameters()); // inits list of "previous" parameters
    setUpMatrices_(sslib);
  }

  virtual Drift* clone() const override
  {
    return new Drift(*this);
  }

  void setUpMatrices_(const SumStatsLibrary& sslib) override;

  void updateMatrices_() override;
};

#endif
