/*
 * Authors: Gustavo V. Barroso
 * Created: 10/08/2022
 * Last modified: 14/06/2023
 *
 */


#include "Mutation.hpp"

// assumes both the infinite sites model as well as equal mutation rates across pops.
void Mutation::setUpMatrices_(const SumStatsLibrary& sslib)
{
  size_t sizeOfBasis = sslib.getSizeOfBasis();
  matrices_.reserve(1);
  std::vector<Eigen::Triplet<double>> coeffs(0);
  coeffs.reserve(sizeOfBasis);

  for(auto it = std::begin(sslib.getBasis()); it != std::end(sslib.getBasis()); ++it)
  {
    int row = it - std::begin(sslib.getBasis());
    int col = -1;

    if((*it)->getPrefix() == "Hl" || (*it)->getPrefix() == "Hr")
    {
      col = sslib.findCompressedIndex(sslib.getDummyIndexUncompressed()); // for a homogeneous system
      coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));
    }

    else if((*it)->getPrefix() == "pi2")
    {
      auto tmp = std::dynamic_pointer_cast<Pi2Moment>(*it);
      assert(tmp != nullptr);

      // introducing 2-locus Het via mutation in right locus (when left already polymorphic)
      col = tmp->getLeftHetStat()->getPosition();
      coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));

      // introducing 2-locus Het via mutation in left locus (when right already polymorphic)
      col = tmp->getRightHetStat()->getPosition();
      coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));
    }

    else if((*it)->getPrefix() != "I" && (*it)->getPrefix() != "DD" && (*it)->getPrefix() != "Dr")
      throw bpp::Exception("Mutation::mis-specified Moment prefix: " + (*it)->getPrefix());
  }

  Eigen::SparseMatrix<double> mat(sizeOfBasis, sizeOfBasis);
  mat.setFromTriplets(std::begin(coeffs), std::end(coeffs));
  mat.makeCompressed();
  mat *= getParameterValue("u");
  matrices_.emplace_back(mat);
  setIdentity_(sizeOfBasis);
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

