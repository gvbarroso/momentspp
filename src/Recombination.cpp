/*
 * Authors: Gustavo V. Barroso
 * Created: 09/08/2022
 * Last modified: 22/08/2022
 *
 */


#include "Recombination.hpp"


void Recombination::setUpMatrices_(size_t matrixSize)
{
  // NOTE for now, this method assumes equal recombination rates across pops.
  // this is why we only access the unique matrix inside matrices_ (using) matrices_[0])
  for(auto it = std::begin(sslib->getStats()); it != std::end(sslib->getStats()); ++it)
  {
    std::string mom = *it->first; // full name of moment
    std::vector<std::string> splitMom = sslib.splitString(mom, "_"); // splits name by underscore

    size_t row = sslib.indexLookup(mom); // recombination matrix only has entries in main diagonal

    if(splitMom[0] == "DD")
      matrices_[0](row, row) -= 2.;

    else if(splitMom[0] == "Dz")
      matrices_[0](row, row) -= 1.;
  }
}

void Recombination::updateMatrices()
{
  std::string paramName = "";

  for(size_t i = 0; i < matrices_.size(); ++i)
  {
    paramName = "r_" + bpp::TexTools::toString(i);

    double prevVal = prevParams_.getParameterValue(paramName);
    double newVal = getParameterValue(paramName); // from within itself

    matrices_[i] *= (newVal / prevVal);
    prevParams_.setParameterValue(paramName, newVal);
  }

  combineMatrices_();
}
