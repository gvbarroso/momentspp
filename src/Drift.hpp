/*
 * Authors: Gustavo V. Barroso
 * Created: 09/08/2022
 * Last modified: 04/09/2025
 *
 */


#ifndef _DRIFT_H_
#define _DRIFT_H_

#include "AbstractOperator.hpp"
#include "SumStatsLibrary.hpp"

class Drift: public AbstractOperator
{

public:
  Drift(const bpp::ParameterList driftParams, const SumStatsLibrary& sslib):
  AbstractOperator(sslib.getPopIndices())
  {
    includeParameters_(driftParams);
    prevParams_.addParameters(getParameters()); // inits list of "previous" parameters
    setUpMatrices_(sslib);
  }

  Drift(const std::vector<long double>& vals, std::shared_ptr<bpp::IntervalConstraint> ic, const SumStatsLibrary& sslib):
  AbstractOperator(sslib.getPopIndices())
  {
    // for each population modeled in the epoch *this operator belongs to, add Ne parameter
    for(size_t i = 0; i < popIndices_.size(); ++i)
      addParameter_(new bpp::Parameter("1/2N_" + bpp::TextTools::toString(popIndices_[i]), vals[i], ic));

    prevParams_.addParameters(getParameters()); // inits list of "previous" parameters
    setUpMatrices_(sslib);
  }

  virtual Drift* clone() const override
  {
    return new Drift(*this);
  }

private:
  int computeDMainDiagContribution_(std::shared_ptr<Moment> mom, size_t id);

  int computeDOffDiagContribution_(std::shared_ptr<Moment> mom, size_t id);

  int computeDrMainDiagContribution_(std::shared_ptr<Moment> mom, size_t id);

  int computeDDMainDiagContribution_(std::shared_ptr<Moment> mom, size_t id);

  int computePi2MainDiagContribution_(std::shared_ptr<Moment> mom, size_t id);

  int computePi2OffDiagContribution_(std::shared_ptr<Moment> mom, size_t id);

  void setUpMatrices_(const SumStatsLibrary& sslib) override;

  void updateMatrices_() override;
};

#endif
