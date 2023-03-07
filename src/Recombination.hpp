/*
 * Authors: Gustavo V. Barroso
 * Created: 09/08/2022
 * Last modified: 07/03/2022
 *
 */


#ifndef _RECOMBINATION_H_
#define _RECOMBINATION_H_

#include "AbstractOperator.hpp"

class Recombination: public AbstractOperator
{

public:
  Recombination(const bpp::ParameterList recParams, const SumStatsLibrary& sslib):
  AbstractOperator(sslib.getNumStats())
  {
    includeParameters_(recParams);
    prevParams_.addParameters(getParameters()); // inits list of "previous" parameters
    setUpMatrices_(sslib);
  }

  Recombination(double initValue, std::shared_ptr<bpp::IntervalConstraint> ic, const SumStatsLibrary& sslib):
  AbstractOperator(sslib.getNumStats())
  {
    //std::shared_ptr<bpp::Parameter> param = std::make_shared<bpp::Parameter>("r_0", initValue, ic);
    //addParameter_(param.get());
    addParameter_(new bpp::Parameter("r_0", initValue, ic));

    prevParams_.addParameters(getParameters()); // inits list of "previous" parameters
    setUpMatrices_(sslib);
  }

  virtual Recombination* clone() const override
  {
    return new Recombination(*this);
  }

  void setUpMatrices_(const SumStatsLibrary& sslib) override;

  void updateMatrices_() override;

};

#endif
