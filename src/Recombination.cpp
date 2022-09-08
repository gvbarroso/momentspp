/*
 * Authors: Gustavo V. Barroso
 * Created: 09/08/2022
 * Last modified: 08/09/2022
 *
 */


#include "Recombination.hpp"


void Recombination::setUpMatrices_(const SumStatsLibrary& sslib)
{
  // for now, this method assumes equal recombination rates across pops.
  matrices_.resize(1);
  eigenDec_.reserve(1);

  matrices_[0] = Eigen::MatrixXd::Zero(sslib.getStats(), sslib.getStats()); // inits to 0 matrix

  for(auto it = std::begin(sslib->getStats()); it != std::end(sslib->getStats()); ++it)
  {
    std::string mom = *it->first; // full name of moment
    std::vector<std::string> splitMom = sslib.splitString(mom, "_"); // splits name by underscore

    size_t row = it - std::begin(sslib->getStats()); // recombination matrix only has entries in main diagonal
    size_t orderD = static_cast<int>(sslib.countInstances(mon, "D"));

    matrices_[i](row, row) = -orderD;
  }

  matrices_[i] *= getParameterValue("r_0");

  EigenDecomposition ed(matrices_[0], exponent_); // is it a problem that mat has zero-columns?
  eigenDec_.emplace_back(ed);
}

void Recombination::updateMatrices_()
{
  std::string paramName = "";

  for(size_t i = 0; i < matrices_.size(); ++i)
  {
    paramName = "r_" + bpp::TexTools::toString(i);

    double prevVal = prevParams_.getParameterValue(paramName);
    double newVal = getParameterValue(paramName);
    double factor = std::pow(newVal / prevVal, exponent_);

    eigenDec_[i].setLambda(eigenDec_[i].lambdaReal() * factor);
  }

  prevParams_.matchParametersValues(getParameters());
}
