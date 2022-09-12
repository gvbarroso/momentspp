/*
 * Authors: Gustavo V. Barroso
 * Created: 10/08/2022
 * Last modified: 12/09/2022
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

  for(auto it = std::begin(sslib.getStats()); it != std::end(sslib.getStats()); ++it)
  {
    std::string mom = *it->first; // full name of focal moment
    std::vector<std::string> splitMom = sslib.splitString(mom, "_"); // splits name by underscore

    std::string p1, p2, p3, p4 = "";

    size_t row = it - std::begin(sslib->getStats()); // row index
    size_t col = 0; // column index

    if(splitMom[0] == "H")
      coefficients.emplace_back(Eigen::Triplet<double>(row, row, 2.)); // main diagonal, introducing one-locus diversity

    else if(splitMom[0] == "pi2")
    {
      std::vector<std::string> splitPops = sslib.splitString(splitMom[1], ";"); // splits name by semi-colon

      p1 = splitPops[0][0]; // i pop
      p2 = splitPops[0][1]; // j pop
      p3 = splitPops[1][0]; // k pop
      p4 = splitPops[1][1]; // l pop

      col = sslib.indexLookup("H_" + p1 + p2);
      coefficients.emplace_back(Eigen::Triplet<double>(row, col, 2.));

      col = sslib.indexLookup("H_" + p3 + p4);
      coefficients.emplace_back(Eigen::Triplet<double>(row, col, 2.));
    }

    else
      coefficients.emplace_back(Eigen::Triplet<double>(row, row, 1.)); // main diagonal, unnaffected terms WARNING 0?
  }

  Eigen::SparseMatrix<double> mat(numStats, numStats);
  mat.setFromTriplets(coefficients);
  mat.makeCompressed();
  mat *= getParameterValue("mu_0");
  matrices_.emplace_back(mat);
}

void Mutation::updateMatrices_()
{
  std::string paramName = "";

  for(size_t i = 0; i < matrices_.size(); ++i)
  {
    paramName = "mu_" + bpp::TextTools::toString(i);

    double prevVal = prevParams_.getParameterValue(paramName);
    double newVal = getParameterValue(paramName);

    matrices_[index] *= (newVal / prevVal);

  }

  prevParams_.matchParametersValues(getParameters());
}

