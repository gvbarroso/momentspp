/*
 * Authors: Gustavo V. Barroso
 * Created: 09/08/2022
 * Last modified: 07/03/2023
 *
 */


#include "Recombination.hpp"

void Recombination::setUpMatrices_(const SumStatsLibrary& sslib)
{
  size_t numStats = sslib.getNumStats();
  // for now, this method assumes equal recombination rates across pops.
  matrices_.reserve(1);
  std::vector<Eigen::Triplet<double>> coeffs(0);
  coeffs.reserve(numStats);

  for(auto it = std::begin(sslib.getMoments()); it != std::end(sslib.getMoments()); ++it)
  {
    int row = it - std::begin(sslib.getMoments()); // recombination matrix only has entries in main diagonal

    if((*it)->getPrefix() == "DD")
      coeffs.emplace_back(Eigen::Triplet<double>(row, row, -2.));

    else if((*it)->getPrefix() == "Dz")
      coeffs.emplace_back(Eigen::Triplet<double>(row, row, -1.));

    else if((*it)->getPrefix() != "I" && (*it)->getPrefix() != "H" && (*it)->getPrefix() != "pi2")
      throw bpp::Exception("Recombination::mis-specified Moment prefix: " + (*it)->getPrefix());
  }

  Eigen::SparseMatrix<double> mat(numStats, numStats);
  mat.setFromTriplets(std::begin(coeffs), std::end(coeffs));
  mat.makeCompressed();
  mat *= getParameterValue("r_0");
  matrices_.emplace_back(mat);
  assembleTransitionMatrix_();
}

void Recombination::updateMatrices_()
{
  double prevVal = prevParams_.getParameterValue("r_0");
  double newVal = getParameterValue("r_0");

  matrices_[0] *= newVal / prevVal;
  assembleTransitionMatrix_();
  prevParams_.matchParametersValues(getParameters());
}
