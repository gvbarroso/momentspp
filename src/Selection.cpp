/*
 * Authors: Gustavo V. Barroso
 * Created: 22/08/2022
 * Last modified: 09/06/2023
 *
 */


#include "Selection.hpp"


void Selection::setUpMatrices_(const SumStatsLibrary& sslib)
{
  size_t sizeOfBasis = sslib.getSizeOfBasis();
  matrices_.reserve(1);
  std::vector<Eigen::Triplet<double>> coeffs(0);
  coeffs.reserve(sizeOfBasis);

  for(auto it = std::begin(sslib.getBasis()); it != std::end(sslib.getBasis()); ++it)
  {
    int row = it - std::begin(sslib.getBasis());
    int col = -1;

    size_t x = (*it)->getFactorPower(); // count of (1-2p) factors on focal moment

    if((*it)->getPrefix() == "DD")
    {
      col = sslib.findCompressedIndex(sslib.findDdIndex((*it)->getPopIndices()[0], (*it)->getPopIndices()[1], x + 1));
      coeffs.emplace_back(Eigen::Triplet<double>(row, col, (2. + x/2.)));

      col = sslib.findCompressedIndex(sslib.findDdIndex((*it)->getPopIndices()[0], (*it)->getPopIndices()[1], x - 1));
      coeffs.emplace_back(Eigen::Triplet<double>(row, col, -x/2.));
    }

    else if((*it)->getPrefix() == "Dr")
    {
      col = sslib.findCompressedIndex(sslib.findDrIndex((*it)->getPopIndices()[0], (*it)->getPopIndices()[1], x + 1));
      coeffs.emplace_back(Eigen::Triplet<double>(row, col, (1. + x/2.)));

      col = sslib.findCompressedIndex(sslib.findDrIndex((*it)->getPopIndices()[0], (*it)->getPopIndices()[1], x - 1));
      coeffs.emplace_back(Eigen::Triplet<double>(row, col, -x/2.));

      col = sslib.findCompressedIndex(sslib.findDdIndex((*it)->getPopIndices()[0], (*it)->getPopIndices()[1], x));
      coeffs.emplace_back(Eigen::Triplet<double>(row, col, -2.));
    }

    else if((*it)->getPrefix() == "Hl")
    {
      col = sslib.findCompressedIndex(sslib.findHetIndex((*it)->getPopIndices()[0], (*it)->getPopIndices()[1], x + 1));
      coeffs.emplace_back(Eigen::Triplet<double>(row, col, (1. + x/2.)));

      col = sslib.findCompressedIndex(sslib.findHetIndex((*it)->getPopIndices()[0], (*it)->getPopIndices()[1], x - 1));
      coeffs.emplace_back(Eigen::Triplet<double>(row, col, -x/2.));
    }

    else if((*it)->getPrefix() == "Hr")
    {
      col = sslib.findCompressedIndex(sslib.findDrIndex((*it)->getPopIndices()[0], (*it)->getPopIndices()[1], 0));
      coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));
    }

    else if((*it)->getPrefix() == "pi2")
    {
      col = sslib.findCompressedIndex(sslib.findPi2Index((*it)->getPopIndices()[0], (*it)->getPopIndices()[1], (*it)->getPopIndices()[2], (*it)->getPopIndices()[3], x + 1));
      coeffs.emplace_back(Eigen::Triplet<double>(row, col, (1. + x/2.)));

      col = sslib.findCompressedIndex(sslib.findPi2Index((*it)->getPopIndices()[0], (*it)->getPopIndices()[1], (*it)->getPopIndices()[2], (*it)->getPopIndices()[3], x - 1));
      coeffs.emplace_back(Eigen::Triplet<double>(row, col, -x/2.));

      col = sslib.findCompressedIndex(sslib.findDrIndex((*it)->getPopIndices()[0], (*it)->getPopIndices()[1], x));
      coeffs.emplace_back(Eigen::Triplet<double>(row, col, 0.25));

      col = sslib.findCompressedIndex(sslib.findDrIndex((*it)->getPopIndices()[0], (*it)->getPopIndices()[1], x + 2));
      coeffs.emplace_back(Eigen::Triplet<double>(row, col, -0.25));
    }

    else if((*it)->getPrefix() != "I" && (*it)->getPrefix() != "DD" && (*it)->getPrefix() != "Dr")
      throw bpp::Exception("Mutation::mis-specified Moment prefix: " + (*it)->getPrefix());
  }

  Eigen::SparseMatrix<double> mat(sizeOfBasis, sizeOfBasis);
  mat.setFromTriplets(std::begin(coeffs), std::end(coeffs));
  mat.makeCompressed();
  mat *= getParameterValue("s");
  matrices_.emplace_back(mat);

  setIdentity_(sizeOfBasis);
  assembleTransitionMatrix_();
}

void Selection::updateMatrices_()
{
  std::string paramName = "s";

  double prevVal = prevParams_.getParameterValue(paramName);
  double newVal = getParameterValue(paramName);

  matrices_.front() *= (newVal / prevVal);
  prevParams_.matchParametersValues(getParameters());
}
