/*
 * Authors: Gustavo V. Barroso
 * Created: 09/08/2022
 * Last modified: 14/06/2023
 *
 */


#include "Recombination.hpp"

// assumes equal recombination rates across pops.
void Recombination::setUpMatrices_(const SumStatsLibrary& sslib)
{
  size_t sizeOfBasis = sslib.getSizeOfBasis();
  matrices_.reserve(1);
  std::vector<Eigen::Triplet<double>> coeffs(0);
  coeffs.reserve(sizeOfBasis);

  for(auto it = std::begin(sslib.getBasis()); it != std::end(sslib.getBasis()); ++it)
  {
    int row = it - std::begin(sslib.getBasis());

    if((*it)->getPrefix() == "DD")
      coeffs.emplace_back(Eigen::Triplet<double>(row, row, -2.));

    else if((*it)->getPrefix() == "Dr")
      coeffs.emplace_back(Eigen::Triplet<double>(row, row, -1.));

    else if((*it)->getPrefix() != "I" && (*it)->getPrefix() != "Hl" && (*it)->getPrefix() != "Hr" && (*it)->getPrefix() != "pi2")
      throw bpp::Exception("Recombination::mis-specified Moment prefix: " + (*it)->getPrefix());
  }

  Eigen::SparseMatrix<double> mat(sizeOfBasis, sizeOfBasis);
  mat.setFromTriplets(std::begin(coeffs), std::end(coeffs));
  mat.makeCompressed();
  mat *= getParameterValue("r");
  matrices_.emplace_back(mat);
  setIdentity_(sizeOfBasis);
  assembleTransitionMatrix_();
}

void Recombination::updateMatrices_()
{
  double prevVal = prevParams_.getParameterValue("r");
  double newVal = getParameterValue("r");

  matrices_[0] *= newVal / prevVal;
  assembleTransitionMatrix_();
  prevParams_.matchParametersValues(getParameters());
}
