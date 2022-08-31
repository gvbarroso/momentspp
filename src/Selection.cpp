/*
 * Authors: Gustavo V. Barroso
 * Created: 22/08/2022
 * Last modified: 31/08/2022
 *
 */


#include "Selection.hpp"


void Selection::setUpMatrices_(const SumStatsLibrary& sslib)
{
  // TODO
}

void Selection::updateMatrices_()
{
  std::string paramName = "";

  for(size_t i = 0; i < matrices_.size(); ++i)
  {
    paramName = "s_" + bpp::TexTools::toString(i);

    double prevVal = prevParams_.getParameterValue(paramName);
    double newVal = getParameterValue(paramName); // from within itself

    solvers_[i].eigenvalues() *= (newVal / prevVal);
    prevParams_.setParameterValue(paramName, newVal);
  }
}
