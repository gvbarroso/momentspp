/*
 * Authors: Gustavo V. Barroso
 * Created: 09/08/2022
 * Last modified: 05/12/2022
 *
 */


#ifndef _DRIFT_H_
#define _DRIFT_H_

#include "AbstractOperator.hpp"
#include "SumStatsLibrary.hpp"

class Drift:
  public AbstractOperator
{

public:
  Drift(const bpp::ParameterList driftParams, const SumStatsLibrary& sslib):
  AbstractOperator(sslib.getNumStats())
  {
    includeParameters_(driftParams);
    prevParams_.addParameters(getParameters()); // inits list of "previous" parameters
    setUpMatrices_(sslib);
  }

  Drift(const std::vector<double>& initValues, std::shared_ptr<bpp::IntervalConstraint> ic, const SumStatsLibrary& sslib):
  AbstractOperator(sslib.getNumStats())
  {
    size_t idx = 0;
    // for each population modeled in the epoch this operator belongs to, add Ne parameter
    for(auto itI = std::begin(sslib.getPopIndices()); itI != std::end(sslib.getPopIndices()); ++itI)
    {
      //std::shared_ptr<bpp::Parameter> param = std::make_shared<bpp::Parameter>("1/2N_" + bpp::TextTools::toString((*itI)), initValue, bpp::Parameter::R_PLUS_STAR);
      //addParameter_(param);
      addParameter_(new bpp::Parameter("1/2N_" + bpp::TextTools::toString(*itI), initValues[idx], ic));
      ++idx;
    }

    prevParams_.addParameters(getParameters()); // inits list of "previous" parameters
    setUpMatrices_(sslib);
  }

  virtual Drift* clone() const override
  {
    return new Drift(*this);
  }

  void setUpMatrices_(const SumStatsLibrary& sslib);

  void updateMatrices_();
};

#endif
