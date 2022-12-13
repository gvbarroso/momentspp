/*
 * Authors: Gustavo V. Barroso
 * Created: 10/08/2022
 * Last modified: 12/12/2022
 *
 */


#include "Mutation.hpp"


void Mutation::setUpMatrices_(const SumStatsLibrary& sslib)
{
  size_t numStats = sslib.getNumStats();

  // for now, this method assumes both the infinite sites model as well as equal mutation rates across pops.
  matrices_.reserve(1);
  std::vector<Eigen::Triplet<double>> coeffs(0);
  coeffs.reserve(numStats);

  for(auto it = std::begin(sslib.getMoments()); it != std::end(sslib.getMoments()); ++it)
  {
    int row = it - std::begin(sslib.getMoments()); // row index
    int col = -1; // inits column index to out-of-bounds

    if(it->getPrefix() == "H") // introducing one-locus diversity
    {
      col = sslib.getDummyMoment().getPosition(); // at column of Dummy Moment
      coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.)); // homogeneous mutation model
    }

    else if(it->getPrefix() == "pi2")
    {
      size_t p1 = it->getPopIndices()[0]; // i pop
      size_t p2 = it->getPopIndices()[1]; // j pop
      size_t p3 = it->getPopIndices()[2]; // k pop
      size_t p4 = it->getPopIndices()[3]; // l pop

      col = sslib.findHetIndex(p1, p2); // introducing 2-locus het via mutation in right locus when left is already polymorphic
      coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));

      col = sslib.findHetIndex(p3, p4); // introducing 2-locus het via mutation in right locus when left is already polymorphic
      coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));
    }
  }

  Eigen::SparseMatrix<double> mat(numStats, numStats);
  mat.setFromTriplets(std::begin(coeffs), std::end(coeffs));
  mat.makeCompressed();
  mat *= getParameterValue("mu_0");
  matrices_.emplace_back(mat);
}

void Mutation::updateMatrices_()
{
  double prevVal = prevParams_.getParameterValue("mu_0");
  double newVal = getParameterValue("mu_0");

  matrices_[0] *= newVal / prevVal;
  prevParams_.matchParametersValues(getParameters());
}

