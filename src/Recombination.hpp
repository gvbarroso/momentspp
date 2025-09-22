/*
 * Authors: Gustavo V. Barroso
 * Created: 09/08/2022
 * Last modified: 18/09/2025
 *
 */

#ifndef _RECOMBINATION_H_
#define _RECOMBINATION_H_

#include "AbstractOperator.hpp"

class Recombination : public AbstractOperator
{

public:
  Recombination(const bpp::ParameterList recParams, const SumStatsLibrary& sslib):
  AbstractOperator(sslib.getPopIndices())
  {
    includeParameters_(recParams);
    prevParams_.addParameters(getParameters()); // inits list of "previous" parameters
    setUpMatrices_(sslib);
  }

  Recombination(const std::vector<double>& initVals, std::shared_ptr<bpp::IntervalConstraint> ic, const SumStatsLibrary& sslib):
  AbstractOperator(sslib.getPopIndices())
  {
    // for each population modeled in the epoch *this operator belongs to, add r parameter
    for(size_t i = 0; i < popIndices_.size(); ++i)
      addParameter_(new bpp::Parameter("r_" + bpp::TextTools::toString(popIndices_[i]), initVals[i], ic));

    prevParams_.addParameters(getParameters()); // inits list of "previous" parameters
    setUpMatrices_(sslib);
  }

  Recombination(const Recombination& other):
  AbstractOperator(other) // invokes base copy constructor
  {
    // no need to call setUpMatrices_ again — matrices_ and transition_ are already cloned
    // parameters are already copied via AbstractOperator's copy constructor
  }

  Recombination* clone() const override
  {
    return new Recombination(*this);
  }

  void setUpMatrices_(const SumStatsLibrary& sslib) override;

  void updateMatrices_() override;
};

#endif
