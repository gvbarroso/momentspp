/*
 * Authors: Gustavo V. Barroso
 * Created: 10/08/2022
 * Last modified: 20/03/2022
 *
 */


#ifndef _MIGRATION_H_
#define _MIGRATION_H_

#include "AbstractOperator.hpp"

class Migration: public AbstractOperator
{

public:
  Migration(const bpp::ParameterList migParams, const SumStatsLibrary& sslib):
  AbstractOperator()
  {
    includeParameters_(migParams);
    prevParams_.addParameters(getParameters()); // inits list of "previous" parameters
    setUpMatrices_(sslib);
  }

  Migration(const std::vector<double>& initValues, std::shared_ptr<bpp::IntervalConstraint> ic, const SumStatsLibrary& sslib):
  AbstractOperator()
  {
    size_t idx = 0;
    for(auto itI = std::begin(sslib.getPopIndices()); itI != std::end(sslib.getPopIndices()); ++itI) // for each population modeled in epoch i
    {
      for(auto itJ = std::begin(sslib.getPopIndices()); itJ != std::end(sslib.getPopIndices()); ++itJ)
      {
        if((*itI) != (*itJ)) // if population indices are different
        {
          //std::shared_ptr<bpp::Parameter> param = std::make_shared<bpp::Parameter>("m_" + bpp::TextTools::toString((*itI)) + bpp::TextTools::toString((*itJ)), initValue, ic);
          //addParameter_(param.get());
          addParameter_(new bpp::Parameter("m_" + bpp::TextTools::toString((*itI)) + "_" + bpp::TextTools::toString((*itJ)), initValues[idx], ic));
          ++idx;
        }
      }
    }

    prevParams_.addParameters(getParameters()); // inits list of "previous" parameters
    setUpMatrices_(sslib);
  }

  virtual Migration* clone() const override
  {
    return new Migration(*this);
  }

  void setUpMatrices_(const SumStatsLibrary& sslib) override;

  void updateMatrices_() override;

  // this is a weird-looking but fun way to get the number of populations P from the raw value of P^2 - P ( == matrices_.size())
  size_t fetchNumPops()
  {
    int numPops = 2; // we want the positive solution of the quadratic equation P^2 - P - matrices_.size() = 0
    int n = static_cast<int>(getParameters().size()); // raw value of P^2 - P

    for(int i = 2; i < n; ++i)
    {
      if(i * (1 - i) == -n)  // guaranteed to find if matrices_.size() was built correctly
      {
        numPops = i;
        break;
      }
    }

    return static_cast<size_t>(numPops);
  }

};

#endif
