/*
 * Authors: Gustavo V. Barroso
 * Created: 09/08/2022
 * Last modified: 15/09/2025
 *
 */

#include <Bpp/Text/TextTools.h>

#include "Recombination.hpp"

void Recombination::setUpMatrices_(const SumStatsLibrary& sslib, bool highPrecision)
{
  const size_t numPops = getParameters().size();
  const size_t basisSize = sslib.getSizeOfBasis();
  matrices_.reserve(numPops);

  const auto& basis = sslib.getBasis();
  const size_t numThreads = omp_get_max_threads();

  for(size_t i = 0; i < numPops; ++i)
  {
    const size_t id = popIndices_[i];
    const std::string paramName = "r_" + bpp::TextTools::toString(id);
    const double recombRate = getParameterValue(paramName);

    // Deduce scalar type
    using Scalar = std::conditional_t<true, mpfr::mpreal, double>;
    if(!highPrecision) using Scalar = double;

    std::vector<std::vector<Eigen::Triplet<Scalar>>> threadTriplets(numThreads);
    std::vector<std::string> invalidPrefixes;

    #pragma omp parallel for
    for(size_t row = 0; row < basisSize; ++row)
    {
      const auto& moment = basis[row];
      const std::string& prefix = moment->getPrefix();
      const size_t tid = omp_get_thread_num();
      auto& localTriplets = threadTriplets[tid];

      if(prefix == "DD")
      {
        size_t f = moment->countInstances(id);
        localTriplets.emplace_back(row, row, Scalar(-f));
      }
      else if(prefix == "Dr" || prefix == "D")
      {
        int f = (moment->getPopIndices()[0] == id);
        localTriplets.emplace_back(row, row, Scalar(-f));
      }
      else if(prefix != "I" && prefix != "Hl" && prefix != "Hr" && prefix != "pi2")
      {
        #pragma omp critical
        invalidPrefixes.push_back(prefix);
      }
    }

    if(!invalidPrefixes.empty())
      throw bpp::Exception("Recombination::mis-specified Moment prefix: " + invalidPrefixes.front());

    // Merge thread-local triplets
    std::vector<Eigen::Triplet<Scalar>> coeffs;
    for(auto& vec : threadTriplets)
      coeffs.insert(coeffs.end(), std::make_move_iterator(vec.begin()), std::make_move_iterator(vec.end()));

    // Create and store matrix
    auto mat = std::make_unique<Matrix<Scalar>>(basisSize, basisSize);
    mat->setFromTriplets(coeffs.begin(), coeffs.end());
    mat->makeCompressed();
    mat->scale(Scalar(recombRate));
    matrices_.emplace_back(std::move(mat));
  }

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
