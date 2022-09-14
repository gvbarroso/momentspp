/*
 * Authors: Gustavo V. Barroso
 * Created: 09/08/2022
 * Last modified: 14/09/2022
 *
 */


#ifndef _RECOMBINATION_H_
#define _RECOMBINATION_H_

#include "Operator.hpp"

class Recombination:
  public Operator
{

public:
  Recombination(std::shared_ptr<bpp::IntervalConstraint> ic, const SumStatsLibrary& ssl):
  Operator()
  {
    addParameter(new bpp::Parameter("r_0", 1e-8, ic));

    prevParams_.addParameters(getParameters()); // inits list of "previous" parameters
    setUpMatrices_(ssl);
  }

  void setUpMatrices_(const SumStatsLibrary& ssl);

  void updateMatrices_();

};

#endif
