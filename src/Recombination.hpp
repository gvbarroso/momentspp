/*
 * Authors: Gustavo V. Barroso
 * Created: 09/08/2022
 * Last modified: 22/09/2022
 *
 */


#ifndef _RECOMBINATION_H_
#define _RECOMBINATION_H_

#include "AbstractOperator.hpp"

class Recombination:
  public AbstractOperator
{

public:
  Recombination(std::shared_ptr<bpp::IntervalConstraint> ic, const SumStatsLibrary& ssl):
  AbstractOperator()
  {
    addParameter_(new bpp::Parameter("r_0", 1e-8, ic));

    prevParams_.addParameters(getParameters()); // inits list of "previous" parameters
    setUpMatrices_(ssl);
  }

  virtual Recombination* clone() const override
  {
    return new Recombination(*this);
  }

  void setUpMatrices_(const SumStatsLibrary& ssl);

  void updateMatrices_();

};

#endif
