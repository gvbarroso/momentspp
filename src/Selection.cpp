/*
 * Authors: Gustavo V. Barroso
 * Created: 22/08/2022
 * Last modified: 16/09/2025
 *
 */

#include "Selection.hpp"

// uses zero-order moment-closure approximation (truncation)
void Selection::setUpMatrices_(const SumStatsLibrary& sslib)
{
  const size_t numPops = getParameters().size();
  const size_t basisSize = sslib.getSizeOfBasis();
  matrices_.reserve(numPops);

  const auto& basis = sslib.getBasis();
  const size_t numThreads = omp_get_max_threads();

  std::visit(
      [&](const auto& mat)
      {
        using Scalar = typename std::decay_t<decltype(mat)>::Scalar;

        for(size_t i = 0; i < numPops; ++i)
        {
          const size_t id = popIndices_[i];
          const std::string paramName = "s_" + bpp::TextTools::toString(id);
          const double s = getParameterValue(paramName);

          std::vector<std::vector<Eigen::Triplet<Scalar>>> threadTriplets(numThreads);

#pragma omp parallel for
          for(size_t row = 0; row < basisSize; ++row)
          {
            int col = -1;

            const auto& moment = basis[row];
            const std::string& prefix = moment->getPrefix();
            const size_t tid = omp_get_thread_num();
            auto& localTriplets = threadTriplets[tid];

            size_t popIdCount = moment->countInstances(
                id); // count of id in moment's name (not counting (1-2p) factors)
            int popIdPower =
                moment->getPopFactorPower(id); // count of (1-2p_x) factors on focal moment

            std::vector<size_t> popIds =
                moment->getPopIndices(); // every moment is guaranteed to have popIds.size() > 1

            if(prefix == "Dr")
            {
              if(popIdCount == 2)
              {
                // Dr contributions
                if(moment->getFactorPower() < sslib.getFactorOrder())
                {
                  std::vector<size_t> factorIds = moment->getFactorIndices();
                  factorIds.push_back(id);

                  col = sslib.findCompressedIndex(sslib.getMoment("Dr", popIds, factorIds));
                  localTriplets.emplace_back(row, col, Scalar(1. + popIdPower / 2.));
                }

                else // truncate
                {
                  std::vector<size_t> factorIds = moment->getFactorIndices();

                  col = sslib.findCompressedIndex(sslib.getMoment("Dr", popIds, factorIds));
                  localTriplets.emplace_back(row, col, Scalar(1. + popIdPower / 2.));
                }

                if(popIdPower > 0)
                {
                  std::vector<size_t> factorIds = moment->getFactorIndices();
                  sslib.dropFactorIds(factorIds, id, 1);

                  col = sslib.findCompressedIndex(sslib.getMoment("Dr", popIds, factorIds));
                  localTriplets.emplace_back(row, col, Scalar(-popIdPower / 2.));
                }

                // DD contributions
                if(moment->getFactorPower() < sslib.getFactorOrder())
                {
                  std::vector<size_t> factorIds = moment->getFactorIndices();

                  col = sslib.findCompressedIndex(sslib.getMoment("DD", popIds, factorIds));
                  localTriplets.emplace_back(row, col, Scalar(-2.));
                }

                else if(popIdPower > 0)
                {
                  std::vector<size_t> factorIds = moment->getFactorIndices();

                  if((moment->getFactorPower() > sslib.getFactorOrder() + 1) && popIdPower > 1)
                    sslib.dropFactorIds(factorIds, id, 2);

                  else
                    sslib.dropFactorIds(factorIds, id, 1);

                  while(factorIds.size() >
                        static_cast<size_t>(sslib.getFactorOrder())) // NOTE truncation
                    factorIds.pop_back();

                  col = sslib.findCompressedIndex(sslib.getMoment("DD", popIds, factorIds));
                  localTriplets.emplace_back(row, col, Scalar(-2.));
                }
              }

              else if(popIdCount == 1)
              {
                if(popIds[0] == id)
                {
                  // +2 because Dr moments include more (1-2px) factors, see
                  // SumStatsLibrary::initMoments_
                  if(moment->getFactorPower() < sslib.getFactorOrder() + 2)
                  {
                    std::vector<size_t> factorIds = moment->getFactorIndices();
                    factorIds.push_back(id);

                    col = sslib.findCompressedIndex(sslib.getMoment("Dr", popIds, factorIds));
                    localTriplets.emplace_back(row, col, Scalar(1. + popIdPower / 2.));
                  }

                  else // truncate
                  {
                    std::vector<size_t> factorIds = moment->getFactorIndices();

                    col = sslib.findCompressedIndex(sslib.getMoment("Dr", popIds, factorIds));
                    localTriplets.emplace_back(row, col, Scalar(1. + popIdPower / 2.));
                  }

                  if(popIdPower > 0)
                  {
                    std::vector<size_t> factorIds = moment->getFactorIndices();
                    sslib.dropFactorIds(factorIds, id, 1);

                    col = sslib.findCompressedIndex(sslib.getMoment("Dr", popIds, factorIds));
                    localTriplets.emplace_back(row, col, Scalar(-popIdPower / 2.));
                  }
                }

                else // if(popIds[1] == id)
                {
                  if(popIdPower > 0)
                  {
                    // +2 because Dr moments include more (1-2px) factors, see
                    // SumStatsLibrary::initMoments_
                    if(moment->getFactorPower() < sslib.getFactorOrder() + 2)
                    {
                      std::vector<size_t> factorIds = moment->getFactorIndices();
                      factorIds.push_back(id);

                      col = sslib.findCompressedIndex(sslib.getMoment("Dr", popIds, factorIds));
                      localTriplets.emplace_back(row, col, Scalar(popIdPower / 2.));
                    }

                    else // truncate
                    {
                      std::vector<size_t> factorIds = moment->getFactorIndices();
                      col = sslib.findCompressedIndex(sslib.getMoment("Dr", popIds, factorIds));
                      localTriplets.emplace_back(row, col, Scalar(popIdPower / 2.));
                    }

                    std::vector<size_t> factorIds = moment->getFactorIndices();
                    sslib.dropFactorIds(factorIds, id, 1);

                    col = sslib.findCompressedIndex(sslib.getMoment("Dr", popIds, factorIds));
                    localTriplets.emplace_back(row, col, Scalar(-popIdPower / 2.));
                  }

                  // DD contributions
                  // +1 because Dr moments include more (1-2px) factors, see
                  // SumStatsLibrary::initMoments_
                  if(moment->getFactorPower() < sslib.getFactorOrder() + 1)
                  {
                    std::vector<size_t> factorIds = moment->getFactorIndices();

                    col = sslib.findCompressedIndex(sslib.getMoment("DD", popIds, factorIds));
                    localTriplets.emplace_back(row, col, Scalar(-2.));
                  }

                  else
                  {
                    if(popIdPower > 0)
                    {
                      std::vector<size_t> factorIds = moment->getFactorIndices();
                      sslib.dropFactorIds(factorIds, id, 1);

                      while(factorIds.size() >
                            static_cast<size_t>(sslib.getFactorOrder())) // NOTE truncation
                        factorIds.pop_back();

                      col = sslib.findCompressedIndex(sslib.getMoment("DD", popIds, factorIds));
                      localTriplets.emplace_back(row, col, Scalar(-2.));
                    }
                  }
                }
              }
            } // ends Dr prefix

#ifdef NAKED_D
            else if(prefix == "D")
            {
              if(moment->getFactorPower() < sslib.getFactorOrder())
              {
                std::vector<size_t> factorIds = moment->getFactorIndices();
                factorIds.push_back(id);

                col = sslib.findCompressedIndex(sslib.getMoment("D", popIds, factorIds));
                localTriplets.emplace_back(row, col, Scalar(1. + popIdPower / 2.));
              }

              else // truncate
              {
                std::vector<size_t> factorIds = moment->getFactorIndices();

                col = sslib.findCompressedIndex(sslib.getMoment("D", popIds, factorIds));
                localTriplets.emplace_back(row, col, Scalar(1. + popIdPower / 2.));
              }

              if(popIdPower > 0)
              {
                std::vector<size_t> factorIds = moment->getFactorIndices();
                sslib.dropFactorIds(factorIds, id, 1);

                col = sslib.findCompressedIndex(sslib.getMoment("D", popIds, factorIds));
                localTriplets.emplace_back(row, col, Scalar(-popIdPower / 2.));
              }
            }
#endif

            else if(prefix == "DD")
            {
              if(moment->getFactorPower() < sslib.getFactorOrder())
              {
                std::vector<size_t> factorIds = moment->getFactorIndices();
                factorIds.push_back(id);

                col = sslib.findCompressedIndex(sslib.getMoment("DD", popIds, factorIds));
                localTriplets.emplace_back(row, col, Scalar(popIdCount + popIdPower / 2.));
              }

              else // truncate
              {
                std::vector<size_t> factorIds = moment->getFactorIndices();

                col = sslib.findCompressedIndex(sslib.getMoment("DD", popIds, factorIds));
                localTriplets.emplace_back(row, col, Scalar(2. + popIdPower / 2.));
              }

              if(popIdPower > 0)
              {
                std::vector<size_t> factorIds = moment->getFactorIndices();
                sslib.dropFactorIds(factorIds, id, 1);

                col = sslib.findCompressedIndex(sslib.getMoment("DD", popIds, factorIds));
                localTriplets.emplace_back(row, col, Scalar(-popIdPower / 2.));
              }
            } // ends DD prefix

            else if(prefix == "Hl")
            {
              if(popIdCount == 2)
              {
                if(moment->getFactorPower() < sslib.getFactorOrder())
                {
                  std::vector<size_t> factorIds = moment->getFactorIndices();
                  factorIds.push_back(id);

                  col = sslib.findCompressedIndex(sslib.getMoment("Hl", popIds, factorIds));
                  localTriplets.emplace_back(row, col, Scalar(1. + popIdPower / 2.));
                }

                else // truncate
                {
                  std::vector<size_t> factorIds = moment->getFactorIndices();

                  col = sslib.findCompressedIndex(sslib.getMoment("Hl", popIds, factorIds));
                  localTriplets.emplace_back(row, col, Scalar(1. + popIdPower / 2.));
                }

                if(popIdPower > 0)
                {
                  std::vector<size_t> factorIds = moment->getFactorIndices();
                  sslib.dropFactorIds(factorIds, id, 1);

                  col = sslib.findCompressedIndex(sslib.getMoment("Hl", popIds, factorIds));
                  localTriplets.emplace_back(row, col, Scalar(-popIdPower / 2.));
                }
              }

              else if(popIdCount == 1)
              {
                localTriplets.emplace_back(row, row, Scalar(std::pow(-1, popIds[0] != id) / 2.));

                if(moment->getFactorPower() < sslib.getFactorOrder())
                {
                  std::vector<size_t> factorIds = moment->getFactorIndices();
                  factorIds.push_back(id);

                  col = sslib.findCompressedIndex(sslib.getMoment("Hl", popIds, factorIds));
                  localTriplets.emplace_back(row, col, Scalar((1. + popIdPower) / 2.));
                }

                else // truncate
                {
                  std::vector<size_t> factorIds = moment->getFactorIndices();

                  col = sslib.findCompressedIndex(sslib.getMoment("Hl", popIds, factorIds));
                  localTriplets.emplace_back(row, col, Scalar((1. + popIdPower) / 2.));
                }

                if(popIdPower > 0)
                {
                  std::vector<size_t> factorIds = moment->getFactorIndices();
                  sslib.dropFactorIds(factorIds, id, 1);

                  col = sslib.findCompressedIndex(sslib.getMoment("Hl", popIds, factorIds));
                  localTriplets.emplace_back(row, col, Scalar(-popIdPower / 2.));
                }
              }

              else if(popIdCount == 0 && popIdPower > 0)
              {
                if(moment->getFactorPower() < sslib.getFactorOrder())
                {
                  std::vector<size_t> factorIds = moment->getFactorIndices();
                  factorIds.push_back(id);

                  col = sslib.findCompressedIndex(sslib.getMoment("Hl", popIds, factorIds));
                  localTriplets.emplace_back(row, col, Scalar(popIdPower / 2.));
                }

                else // truncate
                {
                  std::vector<size_t> factorIds = moment->getFactorIndices();

                  col = sslib.findCompressedIndex(sslib.getMoment("Hl", popIds, factorIds));
                  localTriplets.emplace_back(row, col, Scalar(popIdPower / 2.));
                }

                std::vector<size_t> factorIds = moment->getFactorIndices();
                sslib.dropFactorIds(factorIds, id, 1);

                col = sslib.findCompressedIndex(sslib.getMoment("Hl", popIds, factorIds));
                localTriplets.emplace_back(row, col, Scalar(-popIdPower / 2.));
              }
            } // ends Hl prefix

            else if(prefix == "Hr")
            {
              if(popIdCount == 2)
              {
                std::vector<size_t> factorIds = moment->getFactorIndices();

                col = sslib.findCompressedIndex(sslib.getMoment("Dr", popIds, factorIds));
                localTriplets.emplace_back(row, col, Scalar(1.));
              }

              else if(popIdCount == 1)
              {
                std::vector<size_t> factorIds = moment->getFactorIndices();

                col = sslib.findCompressedIndex(
                    sslib.getMoment("Dr", {id, sslib.fetchOtherId(id)}, factorIds));
                localTriplets.emplace_back(row, col, Scalar(1. / 2.));

                /* cancel out due to right-locus permutation
                #ifdef NAKED_D
                col = sslib.findCompressedIndex(sslib.getMoment("D", { id }, factorIds));
                localTriplets.emplace_back(row, col, 1. / 2.));
                #endif
                */
              }
            } // ends Hr prefix

            else if(prefix == "pi2")
            {
              auto tmpPi2 = std::dynamic_pointer_cast<Pi2Moment>(moment);
              assert(tmpPi2 != nullptr);

              size_t countLeft = tmpPi2->getLeftHetStat()->countInstances(id);
              size_t countRight = tmpPi2->getRightHetStat()->countInstances(id);

              if(countLeft == 1) // has self-contribution, either -s/2 or +s/2
                localTriplets.emplace_back(row, row, Scalar(std::pow(-1, popIds[0] != id) / 2.));

              // pi2 contributions
              if(moment->getFactorPower() < sslib.getFactorOrder())
              {
                std::vector<size_t> factorIds = moment->getFactorIndices();
                factorIds.push_back(id);

                col = sslib.findCompressedIndex(sslib.getMoment("pi2", popIds, factorIds));
                localTriplets.emplace_back(row, col, Scalar((countLeft + popIdPower) / 2.));
              }

              else // NOTE truncation
              {
                std::vector<size_t> factorIds = moment->getFactorIndices();

                col = sslib.findCompressedIndex(sslib.getMoment("pi2", popIds, factorIds));
                localTriplets.emplace_back(row, col, Scalar((countLeft + popIdPower) / 2.));
              }

              if(popIdPower > 0)
              {
                std::vector<size_t> factorIds = moment->getFactorIndices();
                sslib.dropFactorIds(factorIds, id, 1);

                col = sslib.findCompressedIndex(sslib.getMoment("pi2", popIds, factorIds));
                localTriplets.emplace_back(row, col, Scalar(-popIdPower / 2.));
              }

              // contributions from other moments (D and Dr)

              // case: pi2_id_id_id_id
              if((countLeft + countRight) == 4)
              {
                // Dr contributions
                std::vector<size_t> factorIds = moment->getFactorIndices();

                col = sslib.findCompressedIndex(sslib.getMoment("Dr", {id, id}, factorIds));
                localTriplets.emplace_back(row, col, Scalar(1. / 4.));

                if(moment->getFactorPower() < sslib.getFactorOrder())
                {
                  factorIds.push_back(id);
                  factorIds.push_back(id);

                  col = sslib.findCompressedIndex(sslib.getMoment("Dr", {id, id}, factorIds));
                  localTriplets.emplace_back(row, col, Scalar(-1. / 4.));
                }

                else // NOTE truncation
                {
                  factorIds.push_back(id);

                  col = sslib.findCompressedIndex(sslib.getMoment("Dr", {id, id}, factorIds));
                  localTriplets.emplace_back(row, col, Scalar(-1. / 4.));
                }
              }

              // case: pi2_id_id_id_*
              else if(countLeft == 2 && countRight == 1)
              {
                std::vector<size_t> factorIds = moment->getFactorIndices();

                col = sslib.findCompressedIndex(
                    sslib.getMoment("Dr", {popIds[2], popIds[3]}, factorIds));
                localTriplets.emplace_back(row, col, Scalar(1. / 8.));

                if(moment->getFactorPower() < sslib.getFactorOrder())
                {
                  factorIds.push_back(id);
                  factorIds.push_back(id);

                  col = sslib.findCompressedIndex(
                      sslib.getMoment("Dr", {popIds[2], popIds[3]}, factorIds));
                  localTriplets.emplace_back(row, col, Scalar(-1. / 8.));
                }

                else // NOTE truncation
                {
                  factorIds.push_back(id);

                  col = sslib.findCompressedIndex(
                      sslib.getMoment("Dr", {popIds[2], popIds[3]}, factorIds));
                  localTriplets.emplace_back(row, col, Scalar(-1. / 8.));
                }

                /* cancel out due to right-locus permutation
                #ifdef NAKED_D
                factorIds = moment->getFactorIndices(); // reset

                col = sslib.findCompressedIndex(sslib.getMoment("D", { id }, factorIds));
                localTriplets.emplace_back(row, col, Scalar(1. / 4.));

                if(moment->getFactorPower() < sslib.getFactorOrder() - 1)
                {
                  factorIds.push_back(id);
                  factorIds.push_back(id);

                  col = sslib.findCompressedIndex(sslib.getMoment("D", { id }, factorIds));
                  localTriplets.emplace_back(row, col, Scalar(-1. / 4.));
                }

                else if(moment->getFactorPower() < sslib.getFactorOrder())
                {
                  factorIds.push_back(id);

                  col = sslib.findCompressedIndex(sslib.getMoment("D", { id }, factorIds));
                  localTriplets.emplace_back(row, col, Scalar(-1. / 4.));
                }

                else // NOTE heavy truncation makes us collect twice from D_id_{factorIds}
                {
                  col = sslib.findCompressedIndex(sslib.getMoment("D", { id }, factorIds));
                  localTriplets.emplace_back(row, col, Scalar(-1. / 4.));
                }
                #endif
                */
              }

              // case: pi2_id_*_id_id
              else if(popIds[0] == id && countRight == 2)
              {
                std::vector<size_t> factorIds = moment->getFactorIndices();

                col = sslib.findCompressedIndex(sslib.getMoment("Dr", {id, id}, factorIds));
                localTriplets.emplace_back(row, col, Scalar(1. / 4.));

                factorIds.push_back(sslib.fetchOtherId(id));

                col = sslib.findCompressedIndex(sslib.getMoment("Dr", {id, id}, factorIds));
                localTriplets.emplace_back(row, col, Scalar(1. / 4.));

                factorIds.push_back(id);

                col = sslib.findCompressedIndex(sslib.getMoment("Dr", {id, id}, factorIds));
                localTriplets.emplace_back(row, col, Scalar(-1. / 4.));

                sslib.dropFactorIds(factorIds, sslib.fetchOtherId(id), 1);

                col = sslib.findCompressedIndex(sslib.getMoment("Dr", {id, id}, factorIds));
                localTriplets.emplace_back(row, col, Scalar(-1. / 4.));
              }

              // case: pi2_id_*_id_*
              else if(popIds[0] == id && popIds[1] != id && popIds[2] == id && popIds[3] != id)
              {
                std::vector<size_t> factorIds = moment->getFactorIndices();

                col = sslib.findCompressedIndex(sslib.getMoment("Dr", {id, popIds[1]}, factorIds));
                localTriplets.emplace_back(row, col, Scalar(1. / 8.));

                factorIds.push_back(id);

                col = sslib.findCompressedIndex(sslib.getMoment("Dr", {id, popIds[1]}, factorIds));
                localTriplets.emplace_back(row, col, Scalar(-1. / 8.));

                factorIds.push_back(popIds[1]);

                col = sslib.findCompressedIndex(sslib.getMoment("Dr", {id, popIds[1]}, factorIds));
                localTriplets.emplace_back(row, col, Scalar(-1. / 8.));

                sslib.dropFactorIds(factorIds, id, 1);

                col = sslib.findCompressedIndex(sslib.getMoment("Dr", {id, popIds[1]}, factorIds));
                localTriplets.emplace_back(row, col, Scalar(1. / 8.));

#ifdef NAKED_D
                factorIds = moment->getFactorIndices(); // reset

                col = sslib.findCompressedIndex(sslib.getMoment("D", {id}, factorIds));
                localTriplets.emplace_back(row, col, Scalar(1. / 8.));

                factorIds.push_back(id);

                col = sslib.findCompressedIndex(sslib.getMoment("D", {id}, factorIds));
                localTriplets.emplace_back(row, col, Scalar(-1. / 8.));

                factorIds.push_back(popIds[1]);

                col = sslib.findCompressedIndex(sslib.getMoment("D", {id}, factorIds));
                localTriplets.emplace_back(row, col, Scalar(-1. / 8.));

                sslib.dropFactorIds(factorIds, id, 1);

                col = sslib.findCompressedIndex(sslib.getMoment("D", {id}, factorIds));
                localTriplets.emplace_back(row, col, Scalar(1. / 8.));
#endif
              }

              // case: pi2_*_id_id_id
              else if(popIds[0] != id && popIds[1] == id && countRight == 2)
              {
                // Dr contributions
                std::vector<size_t> factorIds = moment->getFactorIndices();

                col = sslib.findCompressedIndex(sslib.getMoment("Dr", {id, id}, factorIds));
                localTriplets.emplace_back(row, col, Scalar(1. / 4.));

                factorIds.push_back(id);

                col = sslib.findCompressedIndex(sslib.getMoment("Dr", {id, id}, factorIds));
                localTriplets.emplace_back(row, col, Scalar(1. / 4.));

                factorIds.push_back(popIds[0]);

                col = sslib.findCompressedIndex(sslib.getMoment("Dr", {id, id}, factorIds));
                localTriplets.emplace_back(row, col, Scalar(-1. / 4.));

                sslib.dropFactorIds(factorIds, id, 1);

                col = sslib.findCompressedIndex(sslib.getMoment("Dr", {id, id}, factorIds));
                localTriplets.emplace_back(row, col, Scalar(-1. / 4.));
              }

              // case: pi2_*_id_id_*
              else if(popIds[0] != id && popIds[1] == id && countRight == 1)
              {
                std::vector<size_t> factorIds = moment->getFactorIndices();

                col = sslib.findCompressedIndex(sslib.getMoment("Dr", {id, popIds[0]}, factorIds));
                localTriplets.emplace_back(row, col, Scalar(1. / 8.));

                factorIds.push_back(id);

                col = sslib.findCompressedIndex(sslib.getMoment("Dr", {id, popIds[0]}, factorIds));
                localTriplets.emplace_back(row, col, Scalar(1. / 8.));

                factorIds.push_back(popIds[0]);

                col = sslib.findCompressedIndex(sslib.getMoment("Dr", {id, popIds[0]}, factorIds));
                localTriplets.emplace_back(row, col, Scalar(-1. / 8.));

                sslib.dropFactorIds(factorIds, id, 1);

                col = sslib.findCompressedIndex(sslib.getMoment("Dr", {id, popIds[0]}, factorIds));
                localTriplets.emplace_back(row, col, Scalar(-1. / 8.));

#ifdef NAKED_D
                factorIds = moment->getFactorIndices(); // reset

                col = sslib.findCompressedIndex(sslib.getMoment("D", {id}, factorIds));
                localTriplets.emplace_back(row, col, Scalar(1. / 8.));

                factorIds.push_back(id);
                col = sslib.findCompressedIndex(sslib.getMoment("D", {id}, factorIds));
                localTriplets.emplace_back(row, col, Scalar(1. / 8.));

                factorIds.push_back(popIds[1]);
                col = sslib.findCompressedIndex(sslib.getMoment("D", {id}, factorIds));
                localTriplets.emplace_back(row, col, Scalar(-1. / 8.));

                sslib.dropFactorIds(factorIds, id, 1);
                col = sslib.findCompressedIndex(sslib.getMoment("D", {id}, factorIds));
                localTriplets.emplace_back(row, col, Scalar(-1. / 8.));
#endif
              }

              // case: pi2_*_*_id_id
              else if(countLeft == 0 && countRight == 2)
              {
                std::vector<size_t> factorIds = moment->getFactorIndices();

                col = sslib.findCompressedIndex(sslib.getMoment("Dr", {id, id}, factorIds));
                localTriplets.emplace_back(row, col, Scalar(1. / 4.));

                factorIds.push_back(popIds[0]);
                factorIds.push_back(popIds[1]);

                col = sslib.findCompressedIndex(sslib.getMoment("Dr", {id, id}, factorIds));
                localTriplets.emplace_back(row, col, Scalar(-1. / 4.));
              }

              // case: pi2_*_*_id_*
              else if(countLeft == 0 && countRight == 1)
              {
                std::vector<size_t> factorIds = moment->getFactorIndices();

                col = sslib.findCompressedIndex(sslib.getMoment("Dr", {id, popIds[0]}, factorIds));
                localTriplets.emplace_back(row, col, Scalar(1. / 8.));

                factorIds.push_back(popIds[0]);
                factorIds.push_back(popIds[1]);

                col = sslib.findCompressedIndex(sslib.getMoment("Dr", {id, popIds[0]}, factorIds));
                localTriplets.emplace_back(row, col, Scalar(-1. / 8.));

#ifdef NAKED_D
                factorIds = moment->getFactorIndices(); // reset

                col = sslib.findCompressedIndex(sslib.getMoment("D", {id}, factorIds));
                localTriplets.emplace_back(row, col, Scalar(1. / 8.));

                factorIds.push_back(popIds[0]);
                factorIds.push_back(popIds[1]);

                col = sslib.findCompressedIndex(sslib.getMoment("D", {id}, factorIds));
                localTriplets.emplace_back(row, col, Scalar(-1. / 8.));
#endif
              }
            } // ends pi2 prefix

            else if(prefix != "I")
              throw bpp::Exception("Selection::mis-specified Moment prefix: " + prefix);
          } // ends parallel loop over basis

          std::vector<Eigen::Triplet<Scalar>> coeffs;
          for(auto& vec : threadTriplets)
            coeffs.insert(coeffs.end(), std::make_move_iterator(vec.begin()), std::make_move_iterator(vec.end()));

          auto mat = std::make_unique<Matrix<Scalar>>(basisSize, basisSize);
          mat->setFromTriplets(coeffs.begin(), coeffs.end());
          mat->makeCompressed();
          mat->scale(Scalar(s));
          matrices_.emplace_back(std::move(mat));
        } // ends loop over pops
      },
      transition_->getMatrixVariant());

  assembleTransitionMatrix_();
}

void Selection::updateMatrices_()
{
  for(size_t i = 0; i < matrices_.size(); ++i)
  {
    size_t id = popIndices_[i];
    std::string paramName = "s_" + bpp::TextTools::toString(id);

    mpfr::mpreal prevVal = prevParams_.getParameterValue(paramName);
    mpfr::mpreal newVal = getParameterValue(paramName);

    if(newVal != prevVal)
      matrices_[i] *= (newVal / prevVal);
  }

  assembleTransitionMatrix_();
  prevParams_.matchParametersValues(getParameters());
}
