/*
 * Authors: Gustavo V. Barroso
 * Created: 10/08/2022
 * Last modified: 03/09/2025
 *
 */


#ifndef _MUTATION_H_
#define _MUTATION_H_

#include "AbstractOperator.hpp"

class Mutation: public AbstractOperator
{

private:
  mpfr::mpreal leftFactor_; // ratio uL / uR

public:
  Mutation(mpfr::mpreal leftFactor, const bpp::ParameterList mutParams, const SumStatsLibrary& sslib):
  AbstractOperator(sslib.getPopIndices()),
  leftFactor_(leftFactor)
  {
    includeParameters_(mutParams);
    prevParams_.addParameters(getParameters()); // inits list of "previous" parameters
    setUpMatrices_(sslib);
  }

  Mutation(mpfr::mpreal leftFactor, const std::vector<mpfr::mpreal>& initVals, std::shared_ptr<bpp::IntervalConstraint> ic, const SumStatsLibrary& sslib):
  AbstractOperator(sslib.getPopIndices()),
  leftFactor_(leftFactor)
  {
    // for each population modeled in the epoch *this operator belongs to, add mu parameter
    for(size_t i = 0; i < popIndices_.size(); ++i)
      addParameter_(new bpp::Parameter("u_" + bpp::TextTools::toString(popIndices_[i]), initVals[i], ic));

    prevParams_.addParameters(getParameters()); // inits list of "previous" parameters
    setUpMatrices_(sslib);
  }

  virtual Mutation* clone() const override
  {
    return new Mutation(*this);
  }

  mpfr::mpreal getLeftFactor()
  {
    return leftFactor_;
  }

  void setUpMatrices_(const SumStatsLibrary& sslib) override;

  void updateMatrices_() override;

};

#endif
