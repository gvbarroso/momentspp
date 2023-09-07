/*
 * Authors: Gustavo V. Barroso
 * Created: 09/08/2022
 * Last modified: 07/09/2023
 *
 */

#include <Bpp/Text/TextTools.h>

#include "Recombination.hpp"

void Recombination::setUpMatrices_(const SumStatsLibrary& sslib)
{
  size_t numPops = getParameters().size();
  size_t sizeOfBasis = sslib.getSizeOfBasis();
  matrices_.reserve(numPops);

  for(size_t i = 0; i < numPops; ++i)
  {
    size_t id = popIndices_[i];
    std::vector<Eigen::Triplet<double>> coeffs(0);
    coeffs.reserve(sizeOfBasis);

    for(auto it = std::begin(sslib.getBasis()); it != std::end(sslib.getBasis()); ++it)
    {
      int row = it - std::begin(sslib.getBasis());

      if((*it)->getPrefix() == "DD")
        coeffs.emplace_back(Eigen::Triplet<double>(row, row, -2.));

      else if((*it)->getPrefix() == "Dr")
        coeffs.emplace_back(Eigen::Triplet<double>(row, row, -1.));

      else if((*it)->getPrefix() == "D")
        coeffs.emplace_back(Eigen::Triplet<double>(row, row, -1.));

      else if((*it)->getPrefix() != "I" && (*it)->getPrefix() != "Hl" && (*it)->getPrefix() != "Hr" && (*it)->getPrefix() != "pi2")
        throw bpp::Exception("Recombination::mis-specified Moment prefix: " + (*it)->getPrefix());
    }

    Eigen::SparseMatrix<double> mat(sizeOfBasis, sizeOfBasis);
    mat.setFromTriplets(std::begin(coeffs), std::end(coeffs));
    mat.makeCompressed();
    mat *= getParameterValue("r_" + bpp::TextTools::toString(id));
    matrices_.emplace_back(mat);
  }

  setIdentity_(sizeOfBasis);
  assembleTransitionMatrix_();
}

void Recombination::updateMatrices_()
{
  for(size_t i = 0; i < matrices_.size(); ++i)
  {
    size_t id = popIndices_[i];
    std::string paramName = "r_" + bpp::TextTools::toString(id);

    double prevVal = prevParams_.getParameterValue(paramName);
    double newVal = getParameterValue(paramName);

    if(newVal != prevVal)
      matrices_[i] *= (newVal / prevVal);
  }

  assembleTransitionMatrix_();
  prevParams_.matchParametersValues(getParameters());
}
