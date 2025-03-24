/*
 * Authors: Gustavo V. Barroso
 * Created: 21/03/2025
 * Last modified: 25/03/2025
 *
 */


#ifndef _NEUTRAL_ADMIXTURE_H_
#define _NEUTRAL_ADMIXTURE_H_

#include "AbstractOperator.hpp"
#include "SumStatsLibrary.hpp"

class NeutralAdmixture: public AbstractOperator
{

private:
  Eigen::MatrixXd littleAdmixMat_; // P x P

public:
  NeutralAdmixture(const bpp::ParameterList admixParams, const SumStatsLibrary& sslib):
  AbstractOperator(sslib.getPopIndices()),
  littleAdmixMat_()
  {
    includeParameters_(admixParams);
    prevParams_.addParameters(getParameters()); // inits list of "previous" parameters
    setUpMatrices_(sslib);
  }

  NeutralAdmixture(const Eigen::MatrixXd& admixMat, const SumStatsLibrary& sslib):
  AbstractOperator(sslib.getPopIndices()),
  littleAdmixMat_(admixMat)
  {
     // unlike model parameters in other operators, f->[0, 1] (ie, no "single event per generation" assumption)
    std::shared_ptr<bpp::IntervalConstraint> ic = std::make_shared<bpp::IntervalConstraint>(0., 1., true, true);
    // for each pair of populations modeled in the epoch to which *this operator belongs
    for(size_t i = 0; i < popIndices_.size(); ++i)
    {
      size_t id = popIndices_[i];

      for(size_t j = 0; j < popIndices_.size(); ++j)
      {
        size_t jd = popIndices_[j];

        if(id != jd && littleAdmixMat_(i, j) > 0.)
        {
          addParameter_(new bpp::Parameter("a_" + bpp::TextTools::toString(id) + "_" + bpp::TextTools::toString(jd), littleAdmixMat_(i, j), ic));
          break; // we don't include g (a.k.a. 1-f) as a model parameter
        }
      }
    }

    prevParams_.addParameters(getParameters());
    setUpMatrices_(sslib);
  }

  virtual NeutralAdmixture* clone() const override
  {
    return new NeutralAdmixture(*this);
  }

  const Eigen::MatrixXd& getLittleAdmixMat()
  {
    return littleAdmixMat_;
  }

private:
  void setUpMatrices_(const SumStatsLibrary& sslib) override;

  void updateMatrices_() override;

  void assembleTransitionMatrix_() override;
};

#endif
