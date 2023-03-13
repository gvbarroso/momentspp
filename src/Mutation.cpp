/*
 * Authors: Gustavo V. Barroso
 * Created: 10/08/2022
 * Last modified: 13/03/2023
 *
 */


#include "Mutation.hpp"

// assumes both the infinite sites model as well as equal mutation rates across pops.
void Mutation::setUpMatrices_(const SumStatsLibrary& sslib)
{
  size_t numStats = sslib.getNumStats();
  matrices_.reserve(1);
  std::vector<Eigen::Triplet<double>> coeffs(0);
  coeffs.reserve(numStats);

  for(auto it = std::begin(sslib.getMoments()); it != std::end(sslib.getMoments()); ++it)
  {
    int row = it - std::begin(sslib.getMoments());
    int col = -1;

    if((*it)->getPrefix() == "H")
    {
      col = sslib.getDummyMoment()->getPosition(); // for a homogeneous system
      coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));
    }

    else if((*it)->getPrefix() == "pi2")
    {
      auto tmp = std::dynamic_pointer_cast<Pi2Moment>(*it);

      if(tmp != nullptr)
      {
        // introducing 2-locus Het via mutation in right locus (when left already polymorphic)
        col = tmp->getLeftHetStat()->getPosition();
        coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));

        // introducing 2-locus Het via mutation in left locus (when right already polymorphic)
        col = tmp->getRightHetStat()->getPosition();
        coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));
      }

      else
        throw bpp::Exception("Mutation::could not downcast pi2 moment: " + (*it)->getName());
    }

    else if((*it)->getPrefix() != "I" && (*it)->getPrefix() != "DD" && (*it)->getPrefix() != "Dz")
      throw bpp::Exception("Mutation::mis-specified Moment prefix: " + (*it)->getPrefix());
  }

  Eigen::SparseMatrix<double> mat(numStats, numStats);
  mat.setFromTriplets(std::begin(coeffs), std::end(coeffs));
  mat.makeCompressed();
  mat *= getParameterValue("u");
  matrices_.emplace_back(mat);
  assembleTransitionMatrix_();
}

void Mutation::updateMatrices_()
{
  double prevVal = prevParams_.getParameterValue("u");
  double newVal = getParameterValue("u");

  matrices_[0] *= newVal / prevVal;
  assembleTransitionMatrix_();
  prevParams_.matchParametersValues(getParameters());
}

