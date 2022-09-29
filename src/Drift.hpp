/*
 * Authors: Gustavo V. Barroso
 * Created: 09/08/2022
 * Last modified: 29/09/2022
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
  Drift(std::shared_ptr<bpp::IntervalConstraint> ic, const SumStatsLibrary& sslib):
  AbstractOperator(sslib.getNumStats())
  {
    double initValue = 1e-4;

    // for each population modeled in the epoch this operator belongs to, add Ne parameter
    for(auto itI = std::begin(sslib.getPopIndices()); itI != std::end(sslib.getPopIndices()); ++itI)
    {
      //std::shared_ptr<bpp::Parameter> param = std::make_shared<bpp::Parameter>("1/N_" + bpp::TextTools::toString((*itI)), initValue, bpp::Parameter::R_PLUS_STAR);
      //addParameter_(param);
      addParameter_(new bpp::Parameter("1/N_" + bpp::TextTools::toString(*itI), initValue, ic));
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
