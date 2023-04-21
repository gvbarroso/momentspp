/*
 * Authors: Gustavo V. Barroso
 * Created: 10/04/2023
 * Last modified: 10/04/2023
 *
 */


#ifndef _ADMIXTURE_H_
#define _ADMIXTURE_H_

#include "AbstractOperator.hpp"
#include "SumStatsLibrary.hpp"

class Admixture: public AbstractOperator
{

public:
  Admixture(const bpp::ParameterList admixParams, const SumStatsLibrary& sslib):
  AbstractOperator(sslib.getNumStats())
  {
    includeParameters_(admixParams);
    prevParams_.addParameters(getParameters()); // inits list of "previous" parameters
    setUpMatrices_(sslib);
  }

  Admixture(const std::vector<double>& initValues, std::shared_ptr<bpp::IntervalConstraint> ic, const SumStatsLibrary& ssl):
  AbstractOperator(sslib.getNumStats())
  {
    size_t idx = 0;
    for(auto itI = std::begin(sslib.getPopIndices()); itI != std::end(sslib.getPopIndices()); ++itI) // for each population modeled in epoch i
    {
      for(auto itJ = std::begin(sslib.getPopIndices()); itJ != std::end(sslib.getPopIndices()); ++itJ)
      {
        if((*itI) != (*itJ)) // if population indices are different
        {
          std::shared_ptr<bpp::Parameter> param = std::make_shared<bpp::Parameter>("m_" + bpp::TextTools::toString((*itI)) + bpp::TextTools::toString((*itJ)), initValues[idx], ic);
          addParameter_(param.get());
          //addParameter_(new bpp::Parameter("f_" + bpp::TextTools::toString((*itI)) + "_" + bpp::TextTools::toString((*itJ)), initValues[idx], ic));
          //++idx;
        }
      }
    }

    prevParams_.addParameters(getParameters()); // inits list of "previous" parameters
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
