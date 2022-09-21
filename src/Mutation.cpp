/*
 * Authors: Gustavo V. Barroso
 * Created: 10/08/2022
 * Last modified: 21/09/2022
 *
 */


#include "Mutation.hpp"


void Mutation::setUpMatrices_(const SumStatsLibrary& sslib)
{
  size_t numStats = sslib.getNumStats();
  // for now, this method assumes both the infinite sites model as well as equal mutation rates across pops.
  matrices_.reserve(1);
  std::vector<Eigen::Triplet<double>> coefficients(0);
  coefficients.reserve(numStats);

  for(auto it = std::begin(sslib.getMoments()); it != std::end(sslib.getMoments()); ++it)
  {
    size_t row = it - std::begin(sslib.getMoments()); // row index
    size_t col = 0; // column index

    if(it->getPrefix() == "H")
      coefficients.emplace_back(Eigen::Triplet<double>(row, row, 2.)); // main diagonal, introducing one-locus diversity

    else if(it->getPrefix() == "pi2")
    {
      size_t p1 = it->getPopIndices()[0]; // i pop
      size_t p2 = it->getPopIndices()[1]; // j pop
      size_t p3 = it->getPopIndices()[2]; // k pop
      size_t p4 = it->getPopIndices()[3]; // l pop

      col = sslib.findHetIndex(p1, p2);
      coefficients.emplace_back(Eigen::Triplet<double>(row, col, 2.));

      col = sslib.findHetIndex(p3, p4);
      coefficients.emplace_back(Eigen::Triplet<double>(row, col, 2.));
    }
  }

  Eigen::SparseMatrix<double> mat(numStats, numStats);
  mat.setFromTriplets(std::begin(coefficients), std::end(coefficients));
  mat.makeCompressed();
  mat *= getParameterValue("mu_0");
  matrices_.emplace_back(mat);
}

void Mutation::updateMatrices_()
{
  for(size_t i = 0; i < matrices_.size(); ++i)
  {
    double prevVal = prevParams_.getParameterValue("mu_0");
    double newVal = getParameterValue("mu_0");

    matrices_[i] *= (newVal / prevVal);
  }

  prevParams_.matchParametersValues(getParameters());
}

