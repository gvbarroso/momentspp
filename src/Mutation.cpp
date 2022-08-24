/*
 * Authors: Gustavo V. Barroso
 * Created: 10/08/2022
 * Last modified: 24/08/2022
 *
 */


#include "Mutation.hpp"


void Mutation::setUpMatrices_(const SumStatsLibrary& sslib)
{
  // NOTE for now, this method assumes both the infinite sites model as well as equal mutation rates across pops.
  matrices_.reserve(1);
  std::vector<Eigen::Triplet<double>> coefficients(0);
  coefficients.reserve(sslib.getNumStats());

  for(auto it = std::begin(sslib->getStats()); it != std::end(sslib->getStats()); ++it)
  {
    std::string mom = *it->first; // full name of focal moment
    std::vector<std::string> splitMom = sslib.splitString(mom, "_"); // splits name by underscore

    std::string p1, p2, p3, p4 = "";

    size_t row = sslib.indexLookup(mom); // row index
    size_t col = 0; // column index

    if(splitMom[0] == "H")
      coefficients.push_back(Eigen::Triplet<double>(row, row, 2.)); // main diagonal, introducing one-locus diversity

    else if(splitMom[0] == "pi2")
    {
      p1 = splitMom[1][0]; // i pop
      p2 = splitMom[1][1]; // j pop
      p3 = splitMom[1][2]; // k pop
      p4 = splitMom[1][3]; // l pop

      col = sslib.indexLookup("H_" + p1 + p2);
      coefficients.push_back(Eigen::Triplet<double>(row, col, 2.));

      col = sslib.indexLookup("H_" + p3 + p4);
      coefficients.push_back(Eigen::Triplet<double>(row, col, 2.));
    }
  }

  Eigen::SparseMatrix<double> mat;
  mat.setFromTriplets(coefficients);
  mat.makeCompressed();
  matrices_.emplace_back(mat);
}

void Mutation::updateMatrices_()
{
  std::string paramName = "";

  for(size_t i = 0; i < matrices_.size(); ++i)
  {
    paramName = "mu_" + bpp::TexTools::toString(i);

    double prevVal = prevParams_.getParameterValue(paramName);
    double newVal = getParameterValue(paramName); // from within itself

    matrices_[i] *= (newVal / prevVal);
  }

  prevParams_.matchParametersValues(getParameters());
  combineMatrices_();
}

