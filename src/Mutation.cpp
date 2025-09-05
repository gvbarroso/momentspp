/*
 * Authors: Gustavo V. Barroso
 * Created: 10/08/2022
 * Last modified: 05/09/2025
 *
 */

#include <Bpp/Text/TextTools.h>

#include "Mutation.hpp"

// assumes both the infinite sites model as well as equal mutation rates across pops.
void Mutation::setUpMatrices_(const SumStatsLibrary& sslib)
{
  const size_t numPops = getParameters().size();
  const size_t sizeOfBasis = sslib.getSizeOfBasis();
  matrices_.reserve(numPops);

  const auto& basis = sslib.getBasis();
  const size_t basisSize = basis.size();

  for (size_t i = 0; i < numPops; ++i)
  {
    const size_t id = popIndices_[i];
    const std::string paramName = "u_" + bpp::TextTools::toString(id);
    const mpfr::mpreal mutationRate = getParameterValue(paramName);

    // prepare thread-local triplet buffers
    const int numThreads = omp_get_max_threads();
    std::vector<std::vector<Eigen::Triplet<mpfr::mpreal>>> threadTriplets(numThreads);

    #pragma omp parallel for
    for (int row = 0; row < static_cast<int>(basisSize); ++row)
    {
      const auto& moment = basis[row];
      const std::string& prefix = moment->getPrefix();
      const int popIdCount = static_cast<int>(moment->countInstances(id));
      const int tid = omp_get_thread_num();
      auto& localTriplets = threadTriplets[tid];

      if(prefix == "Hl" || prefix == "Hr")
      {
        const int col = sslib.findCompressedIndex(sslib.getMoment("I"));
        const mpfr::mpreal factor = (prefix == "Hl") ? leftFactor_ * popIdCount / 2.0 : popIdCount / 2.0;
        localTriplets.emplace_back(row, col, factor);
      }

      else if(prefix == "pi2")
      {
        const auto tmpPi2 = std::dynamic_pointer_cast<Pi2Moment>(moment);
        if (!tmpPi2)
          continue;  // skip invalid cast

        const auto tempLeft = tmpPi2->getLeftHetStat();
        const auto tempRight = tmpPi2->getRightHetStat();

        localTriplets.emplace_back(row, tempLeft->getPosition(), tempLeft->countInstances(id) / 2.0);
        localTriplets.emplace_back(row, tempRight->getPosition(), tempRight->countInstances(id) / 2.0);
      }

      else if(prefix != "I" && prefix != "DD" && prefix != "Dr" && prefix != "D")
      {
        #pragma omp critical
        {
          throw bpp::Exception("Mutation::mis-specified Moment prefix: " + prefix);
        }
      }
    }

    // merge thread-local triplets
    std::vector<Eigen::Triplet<mpfr::mpreal>> coeffs;
    for(auto& vec : threadTriplets)
      coeffs.insert(coeffs.end(), vec.begin(), vec.end());

    Eigen::SparseMatrix<mpfr::mpreal> mat(sizeOfBasis, sizeOfBasis);
    mat.setFromTriplets(coeffs.begin(), coeffs.end());
    mat.makeCompressed();
    mat *= mutationRate;

    matrices_.emplace_back(std::move(mat));
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

    mpfr::mpreal prevVal = prevParams_.getParameterValue(paramName);
    mpfr::mpreal newVal = getParameterValue(paramName);

    if(newVal != prevVal)
      matrices_[i] *= (newVal / prevVal);
  }

  assembleTransitionMatrix_();
  prevParams_.matchParametersValues(getParameters());
}

