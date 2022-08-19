/*
 * Authors: Gustavo V. Barroso
 * Created: 10/08/2022
 * Last modified: 18/08/2022
 *
 */


#include "Mutation.hpp"


void Mutation::setUpMatrix(const SumStatsLibrary& sslib)
{
  // NOTE for now, this method assumes both the infinite sites model as well as equal mutation rates across pops.
  // this is why we only access the unique matrix inside matrices_ (using) matrices_[0])
  for(auto it = std::begin(sslib->getStats()); it != std::end(sslib->getStats()); ++it)
  {
    std::string mom = *it->first; // full name of focal moment
    std::vector<std::string> splitMom = sslib.splitString(mom, "_"); // splits name by underscore

    size_t row = sslib.indexLookup(mom); // row index
    size_t col = 0; // column index

    if(splitMom[0] == "H")
      matrices_[0](row, row) = 2.; // main diagonal, introducing 1-locus diversity

    else if(splitMom[0] == "pi2")
    {
      col = ;
      matrices_[0](row, col) = 2.;
    }

}

void Mutation::update_()
{
  std::string paramName = "";

  for(size_t i = 0; i < matrices_.size(); ++i)
  {
    paramName = "mu_" + bpp::TexTools::toString(i);

    double prevVal = prevParams_.getParameterValue(paramName);
    double newVal = getParameterValue(paramName); // from within itself

    matrices_[i] *= (newVal / prevVal);
    prevParams_.setParameterValue(paramName, newVal);
  }

  combineMatrices_();
}

