/*
 * Authors: Gustavo V. Barroso
 * Created: 10/04/2023
 * Last modified: 03/06/2024
 *
 */


// The Admixture operator on the selection basis (restricted to 2 populations)
#ifndef _ADMIXTURE_H_
#define _ADMIXTURE_H_

#include "AbstractOperator.hpp"
#include "SumStatsLibrary.hpp"

class Admixture: public AbstractOperator
{

private:
  Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> littleAdmixMat_; // 2 x 2

public:
  Admixture(const bpp::ParameterList admixParams, const SumStatsLibrary& sslib):
  AbstractOperator(sslib.getPopIndices()),
  littleAdmixMat_()
  {
    includeParameters_(admixParams);
    prevParams_.addParameters(getParameters()); // inits list of "previous" parameters
    setUpMatrices_(sslib);
  }

  Admixture(const Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic>& admixMat, const SumStatsLibrary& sslib):
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

  virtual Admixture* clone() const override
  {
    return new Admixture(*this);
  }

  const Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic>& getLittleAdmixMat()
  {
    return littleAdmixMat_;
  }

private:
  void setUpMatrices_(const SumStatsLibrary& sslib) override;

  void updateMatrices_() override;

  void assembleTransitionMatrix_() override;
};

#endif
