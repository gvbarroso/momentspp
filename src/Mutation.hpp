/*
 * Authors: Gustavo V. Barroso
 * Created: 10/08/2022
 * Last modified: 24/09/2025
 *
 */

#ifndef _MUTATION_H_
#define _MUTATION_H_

#include "AbstractOperator.hpp"

class Mutation : public AbstractOperator
{

private:
  double leftFactor_; // ratio uL / uR

public:
  Mutation(double leftFactor, const bpp::ParameterList mutParams, const SumStatsLibrary& sslib, bool highPrecision):
  AbstractOperator(sslib.getPopIndices(), highPrecision),
  leftFactor_(leftFactor)
  {
    includeParameters_(mutParams);
    prevParams_.addParameters(getParameters()); // inits list of "previous" parameters
    setUpMatrices_(sslib);
  }

  Mutation(double leftFactor, const std::vector<double>& initVals, std::shared_ptr<bpp::IntervalConstraint> ic, const SumStatsLibrary& sslib, bool highPrecision):
  AbstractOperator(sslib.getPopIndices(), highPrecision),
  leftFactor_(leftFactor)
  {
    // for each population modeled in the epoch *this operator belongs to, add mu parameter
    for(size_t i = 0; i < popIndices_.size(); ++i)
      addParameter_(new bpp::Parameter("u_" + bpp::TextTools::toString(popIndices_[i]), initVals[i], ic));

    prevParams_.addParameters(getParameters()); // inits list of "previous" parameters
    setUpMatrices_(sslib);
  }

  Mutation(const Mutation& other):
  AbstractOperator(other),
  leftFactor_(other.leftFactor_)
  {
    // no need to call setUpMatrices_ again — matrices_ and transition_ are already cloned
    // parameters are already copied via AbstractOperator's copy constructor
  }

  Mutation* clone() const override { return new Mutation(*this); }

  double getLeftFactor()
  {
    return leftFactor_;
  }

  void setUpMatrices_(const SumStatsLibrary& sslib) override;

  void updateMatrices_() override;
};

#endif
