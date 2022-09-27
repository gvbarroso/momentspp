/*
 * Authors: Gustavo V. Barroso
 * Created: 09/08/2022
 * Last modified: 27/09/2022
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
  Drift(const SumStatsLibrary& sslib):
  AbstractOperator(sslib.getNumStats())
  {
    // for each population modeled in the epoch this operator belongs to, add Ne parameter
    for(auto itI = std::begin(sslib.getPopIndices()); itI != std::end(sslib.getPopIndices()); ++itI)
    {
      std::shared_ptr<bpp::Parameter> param = std::make_shared<bpp::Parameter>("N_" + bpp::TextTools::toString((*itI)), 1e+4, bpp::Parameter::R_PLUS_STAR);
      addParameter_(param.get()); // NOTE >= 1e+3?
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
