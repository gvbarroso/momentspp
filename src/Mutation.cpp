/*
 * Authors: Gustavo V. Barroso
 * Created: 10/08/2022
 * Last modified: 26/09/2022
 *
 */


#include "Mutation.hpp"

void Mutation::setUpMatrices_(const SumStatsLibrary& sslib)
{
  size_t numStats = sslib.getNumStats();
  // for now, this method assumes both the infinite sites model as well as equal mutation rates across pops.
  matrices_.reserve(1);
  std::vector<Eigen::Triplet<double>> multCoefficients(0);
  multCoefficients.reserve(numStats);

  // fills in delta matrix w.r.t. two-locus diversity (via matrix multiplication)
  for(auto it = std::begin(sslib.getMoments()); it != std::end(sslib.getMoments()); ++it)
  {
    size_t row = it - std::begin(sslib.getMoments()); // row index
    size_t col = 0; // column index

    if(it->getPrefix() == "pi2")
    {
      size_t p1 = it->getPopIndices()[0]; // i pop
      size_t p2 = it->getPopIndices()[1]; // j pop
      size_t p3 = it->getPopIndices()[2]; // k pop
      size_t p4 = it->getPopIndices()[3]; // l pop

      col = sslib.findHetIndex(p1, p2);
      multCoefficients.emplace_back(Eigen::Triplet<double>(row, col, 2.));

      col = sslib.findHetIndex(p3, p4);
      multCoefficients.emplace_back(Eigen::Triplet<double>(row, col, 2.));
    }
  }

  Eigen::SparseMatrix<double> multMat(numStats, numStats);
  multMat.setFromTriplets(std::begin(multCoefficients), std::end(multCoefficients));
  multMat.makeCompressed();
  multMat *= getParameterValue("mu_0");
  matrices_.emplace_back(multMat);

  std::vector<Eigen::Triplet<double>> addCoefficients(0);
  addCoefficients.reserve(numStats);
  // fills in special delta matrix w.r.t. one-locus diversity (via matrix addition in Epoch::computeExpectedSumStats)
  for(auto it = std::begin(sslib.getMoments()); it != std::end(sslib.getMoments()); ++it)
  {
    size_t row = it - std::begin(sslib.getMoments()); // row index

    if(it->getPrefix() == "H")
      addCoefficients.emplace_back(Eigen::Triplet<double>(row, row, 2.)); // main diagonal, introducing one-locus diversity
  }

  Eigen::SparseMatrix<double> addMat(numStats, numStats);
  addMat.setFromTriplets(std::begin(addCoefficients), std::end(addCoefficients));
  addMat.makeCompressed();
  addMat *= getParameterValue("mu_0");
  oneLocusPi_ = addMat;
}

void Mutation::updateMatrices_()
{
  // single matrix
  double prevVal = prevParams_.getParameterValue("mu_0");
  double newVal = getParameterValue("mu_0");

  matrices_[0] *= (newVal / prevVal);
  oneLocusPi_ *= (newVal / prevVal);

  prevParams_.matchParametersValues(getParameters());
}

