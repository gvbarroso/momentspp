/*
 * Authors: Gustavo V. Barroso
 * Created: 09/08/2022
 * Last modified: 08/09/2025
 *
 */

#include <Bpp/Text/TextTools.h>

#include "Recombination.hpp"

void Recombination::setUpMatrices_(const SumStatsLibrary& sslib, bool highPrecision)
{
  const size_t numPops = getParameters().size();
  const size_t sizeOfBasis = sslib.getSizeOfBasis();
  matrices_.reserve(numPops);

  const auto& basis = sslib.getBasis();
  const size_t basisSize = basis.size();

  for(size_t i = 0; i < numPops; ++i)
  {
    const size_t id = popIndices_[i];
    const std::string paramName = "r_" + bpp::TextTools::toString(id);
    const double recombRate = getParameterValue(paramName);

    // prepare thread-local triplet buffers
    const size_t numThreads = omp_get_max_threads();
    std::vector<std::vector<Eigen::Triplet<double>>> threadTriplets(numThreads);

    #pragma omp parallel for
    for(size_t row = 0; row < static_cast<int>(basisSize); ++row)
    {
      const auto& moment = basis[row];
      const std::string& prefix = moment->getPrefix();
      const size_t tid = omp_get_thread_num();
      auto& localTriplets = threadTriplets[tid];

      if(prefix == "DD")
      {
        size_t f = static_cast<int>(moment->countInstances(id));
        localTriplets.emplace_back(row, row, -f);
      }

      else if(prefix == "Dr" || prefix == "D")
      {
        int f = (moment->getPopIndices()[0] == id);
        localTriplets.emplace_back(row, row, -f);
      }

      else if(prefix != "I" && prefix != "Hl" && prefix != "Hr" && prefix != "pi2")
      {
        #pragma omp critical
        {
          throw bpp::Exception("Recombination::mis-specified Moment prefix: " + prefix);
        }
      }
    }

    // merge thread-local triplets
    std::vector<Eigen::Triplet<double>> coeffs;
    for(auto& vec : threadTriplets)
      coeffs.insert(coeffs.end(), vec.begin(), vec.end());

    if(highPrecision)
    {
      auto mat = std::make_unique<MatrixMPReal>(numStats, numStats);
      mat.setFromTriplets(coeffs.begin(), coeffs.end());
      mat.makeCompressed();
      mat->scale(recombRate);
      matrices_.emplace_back(std::move(mat));
    }

    else
    {
      auto mat = std::make_unique<MatrixDouble>(numStats, numStats);
      mat.setFromTriplets(coeffs.begin(), coeffs.end());
      mat.makeCompressed();
      mat->scale(recombRate);
      matrices_.emplace_back(std::move(mat));
    }
  } // ends loop over populations

  setIdentity_(sizeOfBasis);
  assembleTransitionMatrix_();
}


void Recombination::updateMatrices_()
{
  for(size_t i = 0; i < matrices_.size(); ++i)
  {
    size_t id = popIndices_[i];
    std::string paramName = "r_" + bpp::TextTools::toString(id);

    mpfr::mpreal prevVal = prevParams_.getParameterValue(paramName);
    mpfr::mpreal newVal = getParameterValue(paramName);

    if(newVal != prevVal)
      matrices_[i] *= (newVal / prevVal);
  }

  assembleTransitionMatrix_();
  prevParams_.matchParametersValues(getParameters());
}
