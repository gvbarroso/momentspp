/*
 * Authors: Gustavo V. Barroso
 * Created: 22/08/2022
 * Last modified: 22/08/2022
 *
 */


#include "Selection.hpp"


void Selection::setUpMatrices_(const SumStatsLibrary& sslib)
{
  // TODO
  combineMatrices_();
}

void Selection::updateMatrices()
{
  std::string paramName = "";

  for(size_t i = 0; i < matrices_.size(); ++i)
  {
    paramName = "s_" + bpp::TexTools::toString(i);

    double prevVal = prevParams_.getParameterValue(paramName);
    double newVal = getParameterValue(paramName); // from within itself

    matrices_[i] *= (newVal / prevVal);
    prevParams_.setParameterValue(paramName, newVal);
  }

  combineMatrices_();
}
