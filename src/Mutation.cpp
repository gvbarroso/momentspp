/*
 * Authors: Gustavo V. Barroso
 * Created: 10/08/2022
 * Last modified: 24/09/2025
 *
 */

#include <Bpp/Text/TextTools.h>

#include "Mutation.hpp"

// assumes both the infinite sites model as well as equal mutation rates across pops.
void Mutation::setUpMatrices_(const SumStatsLibrary& sslib)
{
  const size_t numPops = getParameters().size();
  const size_t basisSize = sslib.getSizeOfBasis();
  matrices_.reserve(numPops);

  const auto& basis = sslib.getBasis();
  const size_t numThreads = omp_get_max_threads();

  std::visit(overloaded{[&](auto const& mat)
  {
    using Scalar = typename std::decay_t<decltype(mat)>::Scalar;

    for(size_t i = 0; i < numPops; ++i)
    {
      const size_t id = popIndices_[i];
      const std::string paramName = "u_" + bpp::TextTools::toString(id);
      const double mutationRate = getParameterValue(paramName);

      std::vector<std::vector<Eigen::Triplet<Scalar>>> threadTriplets(numThreads);

      #pragma omp parallel for
      for(size_t row = 0; row < basisSize; ++row)
      {
        const auto& moment = basis[row];
        const std::string& prefix = moment->getPrefix();
        const size_t popIdCount = moment->countInstances(id);
        const size_t tid = omp_get_thread_num();
        auto& localTriplets = threadTriplets[tid];

        if(prefix == "Hl" || prefix == "Hr")
        {
          const size_t col = sslib.findCompressedIndex(sslib.getMoment("I"));
          Scalar factor = (prefix == "Hl") ? Scalar(leftFactor_ * popIdCount / 2.0)
                                            : Scalar(popIdCount / 2.0);
          localTriplets.emplace_back(row, col, factor);
        }

        else if(prefix == "pi2")
        {
          const auto tmpPi2 = std::dynamic_pointer_cast<Pi2Moment>(moment);
          if(!tmpPi2)
          continue;

          const auto tempLeft = tmpPi2->getLeftHetStat();
          const auto tempRight = tmpPi2->getRightHetStat();

          localTriplets.emplace_back(row, tempLeft->getPosition(), Scalar(tempLeft->countInstances(id) / 2.0));
          localTriplets.emplace_back(row, tempRight->getPosition(), Scalar(tempRight->countInstances(id) / 2.0));
        }
      } // ends loop over basis

      std::vector<Eigen::Triplet<Scalar>> coeffs;
      for(auto& vec : threadTriplets)
        coeffs.insert(coeffs.end(), std::make_move_iterator(vec.begin()), std::make_move_iterator(vec.end()));

      auto matrix = std::make_unique<Matrix<Scalar>>(basisSize, basisSize);
      matrix->setFromTriplets(coeffs);
      matrix->makeCompressed();
      matrix->scale(Scalar(mutationRate));

      MatrixEngine::MatrixVariant mv(std::move(*matrix));
      auto engine = std::make_unique<MatrixEngine>(mv, MatrixEngine::VectorVariant{});
      matrices_.emplace_back(std::move(engine));
    } // ends loop over pops
  } // overloaded
  }, transition_->getMatrixVariant());

  std::cout << "Mutation -- assembled individual matrices.\n"
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
      *matrices_[i] *= (newVal / prevVal);
  }

  assembleTransitionMatrix_();
  prevParams_.matchParametersValues(getParameters());
}
