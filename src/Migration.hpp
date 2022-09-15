/*
 * Authors: Gustavo V. Barroso
 * Created: 10/08/2022
 * Last modified: 14/09/2022
 *
 */


#ifndef _MIGRATION_H_
#define _MIGRATION_H_

#include "AbstractOperator.hpp"

class Migration:
  public AbstractOperator
{

public:
  Migration(std::shared_ptr<bpp::IntervalConstraint> ic, const SumStatsLibrary& ssl):
  AbstractOperator()
  {
    double initValue = 1e-8;

    // NOTE the constraint that individual migration rates are "small" guaranteed that the rows
    // of the matrix (m_ij's) sum to 1, with main diagonal entries = 1 - sum of values < 1e=5
    for(auto itI = std::begin(ssl.getPopMap()); itI != std::end(ssl.getPopMap()); ++itI) // for each population modeled in epoch i
    {
      for(auto itJ = std::begin(ssl.getPopMap()); itJ != std::end(ssl.getPopMap()); ++itJ)
      {
        if((*itI).first != (*itJ).first) // if population indices are different
          migPl.addParameter(new bpp::Parameter("m_" + bpp::TextTools::toString((*itI).first) + bpp::TextTools::toString((*itJ).first), initValue, ic));
      }
    }

    prevParams_.addParameters(getParameters()); // inits list of "previous" parameters
    setUpMatrices_(ssl);
  }

  void setUpMatrices_(const SumStatsLibrary& ssl);

  void updateMatrices_();

  // this is a weird-looking but fun way to get the number of populations P from the raw value of P^2 - P ( == matrices_.size())
  size_t fetchNumPops()
  {
    int numPops = 0; // we want the positive solution of the quadratic equation P^2 - P - matrices_.size() = 0
    int n = static_cast<int>(matrices_.size()); // raw value of P^2 - P

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
