/*
 * Authors: Gustavo V. Barroso
 * Created: 09/08/2022
 * Last modified: 07/09/2022
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

  Eigen::EigenSolver es(matrices_[0]); // is it a problem that mat has zero-columns?
  eigenDec_.emplace_back(es);
}

void Recombination::updateMatrices_()
{
  std::string paramName = "";

  for(size_t i = 0; i < matrices_.size(); ++i)
  {
    paramName = "r_" + bpp::TexTools::toString(i);

    double prevVal = prevParams_.getParameterValue(paramName);
    double newVal = getParameterValue(paramName); // from within itself

    eigenDec_[i].setLambda(eigenDec_[i].lambda() * (newVal / prevVal));
  }

  prevParams_.matchParametersValues(getParameters());
}
