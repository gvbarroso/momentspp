/*
 * Authors: Gustavo V. Barroso
 * Created: 10/08/2022
 * Last modified: 06/06/2024
 *
 */

#include <Bpp/Text/TextTools.h>

#include "Mutation.hpp"

// assumes both the infinite sites model as well as equal mutation rates across pops.
void Mutation::setUpMatrices_(const SumStatsLibrary& sslib)
{
  size_t numPops = getParameters().size();
  size_t sizeOfBasis = sslib.getSizeOfBasis();
  matrices_.reserve(numPops);

  for(size_t i = 0; i < numPops; ++i)
  {
    size_t id = popIndices_[i];
    std::vector<Eigen::Triplet<long double>> coeffs(0);
    coeffs.reserve(sizeOfBasis);

    for(auto it = std::begin(sslib.getBasis()); it != std::end(sslib.getBasis()); ++it)
    {
      int row = it - std::begin(sslib.getBasis());
      int col = -1;

      int popIdCount = static_cast<int>((*it)->countInstances(id));

      if((*it)->getPrefix() == "Hl")
      {
        col = sslib.findCompressedIndex(sslib.getMoment("I")); // for a homogeneous system
        coeffs.emplace_back(Eigen::Triplet<long double>(row, col, leftFactor_ * popIdCount / 2.));
      }

      else if((*it)->getPrefix() == "Hr")
      {
        col = sslib.findCompressedIndex(sslib.getMoment("I")); // for a homogeneous system
        coeffs.emplace_back(Eigen::Triplet<long double>(row, col, popIdCount / 2.)); // ie, rightFactor_ == 1
      }

      else if((*it)->getPrefix() == "pi2")
      {
        auto tmpPi2 = std::dynamic_pointer_cast<Pi2Moment>(*it);
        assert(tmpPi2 != nullptr);

        // introducing 2-locus Het via mutation in right locus (when left already polymorphic)
        auto tempLeft = tmpPi2->getLeftHetStat();
        size_t leftMult = tempLeft->countInstances(id);
        col = tempLeft->getPosition();
        coeffs.emplace_back(Eigen::Triplet<long double>(row, col, leftMult / 2.));

        // introducing 2-locus Het via mutation in left locus (when right already polymorphic)
        auto tempRight = tmpPi2->getRightHetStat();
        size_t rightMult = tempRight->countInstances(id);
        col = tempRight->getPosition();
        coeffs.emplace_back(Eigen::Triplet<long double>(row, col, rightMult / 2.));
      }

      else if((*it)->getPrefix() != "I" && (*it)->getPrefix() != "DD" && (*it)->getPrefix() != "Dr" && (*it)->getPrefix() != "D")
        throw bpp::Exception("Mutation::mis-specified Moment prefix: " + (*it)->getPrefix());
    }

    Eigen::SparseMatrix<long double> mat(sizeOfBasis, sizeOfBasis);
    mat.setFromTriplets(std::begin(coeffs), std::end(coeffs));
    mat.makeCompressed();
    mat *= getParameterValue("u_" + bpp::TextTools::toString(id));
    matrices_.emplace_back(mat);
  }

  setIdentity_(sizeOfBasis);
  assembleTransitionMatrix_();
}

void Mutation::updateMatrices_()
{
  for(size_t i = 0; i < matrices_.size(); ++i)
  {
    size_t id = popIndices_[i];
    std::string paramName = "u_" + bpp::TextTools::toString(id);

    long double prevVal = prevParams_.getParameterValue(paramName);
    long double newVal = getParameterValue(paramName);

    if(newVal != prevVal)
      matrices_[i] *= (newVal / prevVal);
  }

  assembleTransitionMatrix_();
  prevParams_.matchParametersValues(getParameters());
}

