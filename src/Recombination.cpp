/*
 * Authors: Gustavo V. Barroso
 * Created: 09/08/2022
 * Last modified: 24/08/2022
 *
 */


#include "Recombination.hpp"


void Recombination::setUpMatrices_(const SumStatsLibrary& sslib)
{
  // for now, this method assumes equal recombination rates across pops.
  matrices_.reserve(1);
  std::vector<Eigen::Triplet<double>> coefficients(0);
  coefficients.reserve(sslib.getNumStats());

  for(auto it = std::begin(sslib->getStats()); it != std::end(sslib->getStats()); ++it)
  {
    std::string mom = *it->first; // full name of moment
    std::vector<std::string> splitMom = sslib.splitString(mom, "_"); // splits name by underscore

    size_t row = it - std::begin(sslib->getStats()); // recombination matrix only has entries in main diagonal
    size_t orderD = static_cast<int>(sslib.countInstances(mon, "D"));

    coefficients.push_back(Eigen::Triplet<double>(row, row, -orderD));
  }

  Eigen::SparseMatrix<double> mat;
  mat.setFromTriplets(coefficients);
  mat.makeCompressed();
  matrices_.emplace_back(mat);

  Eigen::EigenSolver es(mat); // is it a problem that mat has zero-columns?
  solvers_.emplace_back(es); // TODO check if it's more efficient to store eigenvectors, eigenvalues, and eigenvectors ^ (-1)
}

void Recombination::updateMatrices_()
{
  std::string paramName = "";

  for(size_t i = 0; i < matrices_.size(); ++i)
  {
    paramName = "r_" + bpp::TexTools::toString(i);

    double prevVal = prevParams_.getParameterValue(paramName);
    double newVal = getParameterValue(paramName); // from within itself

    solvers_[i].eigenvalues() *= (newVal / prevVal);
  }

  prevParams_.matchParametersValues(getParameters());
}
