/*
 * Authors: Gustavo V. Barroso
 * Created: 10/08/2022
 * Last modified: 12/12/2022
 *
 */


#ifndef _MUTATION_H_
#define _MUTATION_H_

#include "AbstractOperator.hpp"

class Mutation:
  public AbstractOperator
{

public:
  Mutation(const bpp::ParameterList mutParams, const SumStatsLibrary& sslib):
  AbstractOperator(sslib.getNumStats())
  {
    includeParameters_(mutParams);
    prevParams_.addParameters(getParameters()); // inits list of "previous" parameters
    setUpMatrices_(sslib);
  }

  Mutation(double initValue, std::shared_ptr<bpp::IntervalConstraint> ic, const SumStatsLibrary& sslib):
  AbstractOperator(sslib.getNumStats())
  {
    //std::shared_ptr<bpp::Parameter> param = std::make_shared<bpp::Parameter>("u", initValue, ic);
    //addParameter_(param.get());
    addParameter_(new bpp::Parameter("u", initValue, ic));

    prevParams_.addParameters(getParameters()); // inits list of "previous" parameters
    setUpMatrices_(sslib);
  }

  virtual Mutation* clone() const override
  {
    return new Mutation(*this);
  }

  void setUpMatrices_(const SumStatsLibrary& sslib) override;

  void updateMatrices_() override;

};

#endif
