/*
 * Authors: Gustavo V. Barroso
 * Created: 09/08/2022
 * Last modified: 23/08/2022
 *
 */


#ifndef _RECOMBINATION_H_
#define _RECOMBINATION_H_

#include "Operator.hpp"

class Recombination:
  public Operator
{

public:
  Recombination(const SumStatsLibrary& ssl):
  Operator()
  {
    bpp::addParameter_(new bpp::Parameter("r_0", 1., Parameter::R_PLUS_STAR)); // NOTE could be 1 r' per genomic window etc

    prevParams_.addParameters(getParameters()); // inits list of "previous" parameters
    setUpMatrices_(ssl);

    setUpMatrices_(ssl);
  }

  Recombination(const bpp::ParameterList& params, const SumStatsLibrary& ssl):
  Operator()
  {
    bpp::addParameter_(new bpp::Parameter("r_0", 1., Parameter::R_PLUS_STAR));

    matchParametersValues(params);
    prevParams_.addParameters(getParameters()); // inits list of "previous" parameters
    setUpMatrices_(ssl);
  }

};

#endif
