/*
 * Authors: Gustavo V. Barroso
 * Created: 10/04/2023
 * Last modified: 24/04/2023
 *
 */


#ifndef _ADMIXTURE_H_
#define _ADMIXTURE_H_

#include "AbstractOperator.hpp"
#include "SumStatsLibrary.hpp"

class Admixture: public AbstractOperator
{

private:
  Eigen::MatrixXd littleAdmixMat_; // P x P

public:
  Admixture(const bpp::ParameterList admixParams, const SumStatsLibrary& sslib):
  AbstractOperator(sslib.getNumStats()),
  littleAdmixMat_()
  {
    includeParameters_(admixParams);
    prevParams_.addParameters(getParameters()); // inits list of "previous" parameters
    setUpMatrices_(sslib);
  }

  Admixture(const Eigen::MatrixXd& admixMat, const SumStatsLibrary& sslib):
  AbstractOperator(sslib.getPopIndices()),
  littleAdmixMat_(admixMat)
  {
    // for each pair of populations modeled in the epoch to which *this operator belongs
    for(size_t i = 0; i < popIndices_.size(); ++i)
    {
      size_t id = popIndices_[i];

      for(size_t j = 0; j < popIndices_.size(); ++j)
      {
        size_t jd = popIndices_[j];

        if(id != jd)
          addParameter_(new bpp::Parameter("a_" + bpp::TextTools::toString(id) + "_" + bpp::TextTools::toString(jd), littleAdmixMat_(i, j)));
      }
    }

    prevParams_.addParameters(getParameters());
    setUpMatrices_(sslib);
  }

  virtual Admixture* clone() const override
  {
    return new Admixture(*this);
  }

private:
  void setUpMatrices_(const SumStatsLibrary& sslib) override;

  void updateMatrices_() override;

};

#endif
