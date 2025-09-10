/*
 * Authors: Gustavo V. Barroso
 * Created: 10/08/2022
 * Last modified: 08/09/2025
 *
 */

#include <Bpp/Text/TextTools.h>

#include "Mutation.hpp"

// assumes both the infinite sites model as well as equal mutation rates across pops.
void Mutation::setUpMatrices_(const SumStatsLibrary& sslib, bool highPrecision)
{
  const size_t numPops = getParameters().size();
  const size_t sizeOfBasis = sslib.getSizeOfBasis();
  matrices_.reserve(numPops);

  const auto& basis = sslib.getBasis();
  const size_t basisSize = basis.size();

  for(size_t i = 0; i < numPops; ++i)
  {
    const size_t id = popIndices_[i];
    const std::string paramName = "u_" + bpp::TextTools::toString(id);
    const double mutationRate = getParameterValue(paramName);

    // prepare thread-local triplet buffers
    const size_t numThreads = omp_get_max_threads();
    std::vector<std::vector<Eigen::Triplet<double>>> threadTriplets(numThreads);

    #pragma omp parallel for
    for(size_t row = 0; row < static_cast<size_t>(basisSize); ++row)
    {
      const auto& moment = basis[row];
      const std::string& prefix = moment->getPrefix();
      const size_t popIdCount = static_cast<size_t>(moment->countInstances(id));
      const size_t tid = omp_get_thread_num();
      auto& localTriplets = threadTriplets[tid];

      if(prefix == "Hl" || prefix == "Hr")
      {
        const size_t col = sslib.findCompressedIndex(sslib.getMoment("I"));
        const double factor = (prefix == "Hl") ? leftFactor_ * popIdCount / 2.0 : popIdCount / 2.0;
        localTriplets.emplace_back(row, col, factor);
      }

      else if(prefix == "pi2")
      {
        const auto tmpPi2 = std::dynamic_pointer_cast<Pi2Moment>(moment);
        if(!tmpPi2)
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
    std::vector<Eigen::Triplet<double>> coeffs;
    for(auto& vec : threadTriplets)
      coeffs.insert(coeffs.end(), vec.begin(), vec.end());

    if(highPrecision)
    {
      auto mat = std::make_unique<MatrixMPReal>(numStats, numStats);
      mat.setFromTriplets(coeffs.begin(), coeffs.end());
      mat.makeCompressed();
      mat->scale(mutationRate);
      matrices_.emplace_back(std::move(mat));
    }

    else
    {
      auto mat = std::make_unique<MatrixDouble>(numStats, numStats);
      mat.setFromTriplets(coeffs.begin(), coeffs.end());
      mat.makeCompressed();
      mat->scale(mutationRate);
      matrices_.emplace_back(std::move(mat));
    }
  } // ends loop over populations

  setIdentity_(sizeOfBasis);
  assembleTransitionMatrix_();
}

void Mutation::updateMatrices_()
{
  for(size_t i = 0; i < matrices_.size(); ++i)
  {
    size_t id = popIndices_[i];
    std::string paramName = "u_" + bpp::TextTools::toString(id);

    double prevVal = prevParams_.getParameterValue(paramName);
    double newVal = getParameterValue(paramName);

    if(newVal != prevVal)
      matrices_[i]->scale(newVal / prevVal);
  }

  assembleTransitionMatrix_();
  prevParams_.matchParametersValues(getParameters());
}

