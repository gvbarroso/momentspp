/*
 * Authors: Gustavo V. Barroso
 * Created: 10/08/2022
 * Last modified: 07/09/2022
 *
 */


#include "Mutation.hpp"


void Mutation::setUpMatrices_(const SumStatsLibrary& sslib)
{
  // for now, this method assumes both the infinite sites model as well as equal mutation rates across pops.
  matrices_.resize(1);
  solvers_.reserve(1);

  matrices_[0] = Eigen::MatrixXd::Zero(sslib.getStats(), sslib.getStats()); // inits to 0 matrix

  for(auto it = std::begin(sslib->getStats()); it != std::end(sslib->getStats()); ++it)
  {
    std::string mom = *it->first; // full name of focal moment
    std::vector<std::string> splitMom = sslib.splitString(mom, "_"); // splits name by underscore

    std::string p1, p2, p3, p4 = "";

    size_t row = it - std::begin(sslib->getStats()); // row index
    size_t col = 0; // column index

    if(splitMom[0] == "H")
      matrices_[0](row, row) = 2.; // main diagonal, introducing one-locus diversity

    else if(splitMom[0] == "pi2")
    {
      std::vector<std::string> splitPops = sslib.splitString(splitMom[1], ";"); // splits name by semi-colon

      p1 = splitPops[0][0]; // i pop
      p2 = splitPops[0][1]; // j pop
      p3 = splitPops[1][0]; // k pop
      p4 = splitPops[1][1]; // l pop

      col = sslib.indexLookup("H_" + p1 + p2);
      matrices_[0](row, col) = 2.;

      col = sslib.indexLookup("H_" + p3 + p4);
      matrices_[0](row, col) = 2.;
    }

    else
      matrices_[0](row, row) = 1.; // main diagonal, unnaffected terms
  }

  Eigen::EigenSolver es(matrices_[0]); // is it a problem that mat has zero-columns?
  eigenDec_.emplace_back(es);
}

void Mutation::updateMatrices_()
{
  std::string paramName = "";

  for(size_t i = 0; i < matrices_.size(); ++i)
  {
    paramName = "mu_" + bpp::TexTools::toString(i);

    double prevVal = prevParams_.getParameterValue(paramName);
    double newVal = getParameterValue(paramName); // from within itself

    eigenDec_[i].setLambda(eigenDec_[i].lambda() * (newVal / prevVal));
  }

  prevParams_.matchParametersValues(getParameters());
}

