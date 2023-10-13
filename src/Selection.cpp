/*
 * Authors: Gustavo V. Barroso
 * Created: 22/08/2022
 * Last modified: 13/10/2023
 *
 */


#include "Selection.hpp"


// uses zero-order moment-closure approximation (truncation)
void Selection::setUpMatrices_(const SumStatsLibrary& sslib)
{
  size_t numPops = getParameters().size();
  size_t sizeOfBasis = sslib.getSizeOfBasis();
  matrices_.reserve(numPops);

  for(size_t i = 0; i < numPops; ++i)
  {
    size_t id = popIndices_[i];
    std::vector<Eigen::Triplet<double>> coeffs(0);
    coeffs.reserve(sizeOfBasis);

    for(auto it = std::begin(sslib.getBasis()); it != std::end(sslib.getBasis()); ++it)
    {
      int row = it - std::begin(sslib.getBasis());
      int col = -1;

      size_t popIdCount = (*it)->countInstances(id); // count of id in moment's name (not counting (1-2p) factors)
      int popIdPower = (*it)->getPopFactorPower(id); // count of (1-2p_x) factors on focal moment

      std::vector<size_t> popIds = (*it)->getPopIndices();

      if((*it)->getPrefix() == "Dr")
      {
        if(popIdCount == 2)
        {
          // Dr contributions
          if((*it)->getFactorPower() < sslib.getFactorOrder())
          {
            std::vector<size_t> factorIds = (*it)->getFactorIndices();
            factorIds.push_back(id);

            col = sslib.findCompressedIndex(sslib.getMoment("Dr", popIds, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, (1. + popIdPower / 2.)));
          }

          else // truncate
          {
            std::vector<size_t> factorIds = (*it)->getFactorIndices();

            col = sslib.findCompressedIndex(sslib.getMoment("Dr", popIds, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, (1. + popIdPower / 2.)));
          }

          if(popIdPower > 0)
          {
            std::vector<size_t> factorIds = (*it)->getFactorIndices();
            sslib.dropFactorIds(factorIds, id, 1);

            col = sslib.findCompressedIndex(sslib.getMoment("Dr", popIds, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, -popIdPower / 2.));
          }

          // DD contributions
          if((*it)->getFactorPower() < sslib.getFactorOrder())
          {
            std::vector<size_t> factorIds = (*it)->getFactorIndices();

            col = sslib.findCompressedIndex(sslib.getMoment("DD", popIds, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, -2.));
          }

          else if(popIdPower > 0)
          {
            std::vector<size_t> factorIds = (*it)->getFactorIndices();

            if(((*it)->getFactorPower() > sslib.getFactorOrder() + 1) && popIdPower > 1)
              sslib.dropFactorIds(factorIds, id, 2);

            else
              sslib.dropFactorIds(factorIds, id, 1);

            while(factorIds.size() > static_cast<size_t>(sslib.getFactorOrder())) // NOTE truncation
              factorIds.pop_back();

            col = sslib.findCompressedIndex(sslib.getMoment("DD", popIds, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, -2.));
          }
        }

        else if(popIdCount == 1)
        {
          if(popIds[0] == id)
          {
            if((*it)->getFactorPower() < sslib.getFactorOrder())
            {
              std::vector<size_t> factorIds = (*it)->getFactorIndices();
              factorIds.push_back(id);

              col = sslib.findCompressedIndex(sslib.getMoment("Dr", popIds, factorIds));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, (1. + popIdPower / 2.)));
            }

            else // truncate
            {
              std::vector<size_t> factorIds = (*it)->getFactorIndices();

              col = sslib.findCompressedIndex(sslib.getMoment("Dr", popIds, factorIds));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, (1. + popIdPower / 2.)));
            }

            if(popIdPower > 0)
            {
              std::vector<size_t> factorIds = (*it)->getFactorIndices();
              sslib.dropFactorIds(factorIds, id, 1);

              col = sslib.findCompressedIndex(sslib.getMoment("Dr", popIds, factorIds));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, -popIdPower / 2.));
            }
          }

          else // if(popIds[1] == id)
          {
            if((*it)->getFactorPower() < sslib.getFactorOrder())
            {
              std::vector<size_t> factorIds = (*it)->getFactorIndices();
              factorIds.push_back(id);

              col = sslib.findCompressedIndex(sslib.getMoment("Dr", popIds, factorIds));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, (1. + popIdPower / 2.)));
            }

            else // truncate
            {
              std::vector<size_t> factorIds = (*it)->getFactorIndices();
              col = sslib.findCompressedIndex(sslib.getMoment("Dr", popIds, factorIds));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, (1. + popIdPower / 2.)));
            }

            if(popIdPower > 0)
            {
              std::vector<size_t> factorIds = (*it)->getFactorIndices();
              sslib.dropFactorIds(factorIds, id, 1);

              col = sslib.findCompressedIndex(sslib.getMoment("Dr", popIds, factorIds));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, -popIdPower / 2.));
            }

            // DD contributions
            if((*it)->getFactorPower() < sslib.getFactorOrder())
            {
              std::vector<size_t> factorIds = (*it)->getFactorIndices();

              col = sslib.findCompressedIndex(sslib.getMoment("DD", popIds, factorIds));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, -2.));
            }

            else if(popIdPower > 0)
            {
              std::vector<size_t> factorIds = (*it)->getFactorIndices();
              sslib.dropFactorIds(factorIds, id, 1);

              while(factorIds.size() > static_cast<size_t>(sslib.getFactorOrder())) // NOTE truncation
                factorIds.pop_back();

              col = sslib.findCompressedIndex(sslib.getMoment("DD", popIds, factorIds));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, -2.));
            }
          }
        }
      }

      #ifdef NAKED_D
      else if((*it)->getPrefix() == "D")
      {
        if((*it)->getFactorPower() < sslib.getFactorOrder())
        {
          std::vector<size_t> factorIds = (*it)->getFactorIndices();
          factorIds.push_back(id);

          col = sslib.findCompressedIndex(sslib.getMoment("D", popIds, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, (1. + popIdPower/2.)));
        }

        else // truncate
        {
          std::vector<size_t> factorIds = (*it)->getFactorIndices();

          col = sslib.findCompressedIndex(sslib.getMoment("D", popIds, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, (1. + popIdPower/2.)));
        }

        if(popIdPower > 0)
        {
          std::vector<size_t> factorIds = (*it)->getFactorIndices();
          sslib.dropFactorIds(factorIds, id, 1);

          col = sslib.findCompressedIndex(sslib.getMoment("D", popIds, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, -popIdPower/2.));
        }
      }
      #endif

      else if((*it)->getPrefix() == "DD")
      {
        if((*it)->getFactorPower() < sslib.getFactorOrder())
        {
          std::vector<size_t> factorIds = (*it)->getFactorIndices();
          factorIds.push_back(id);

          col = sslib.findCompressedIndex(sslib.getMoment("DD", popIds, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, (popIdCount + popIdPower / 2.)));
        }

        else // truncate
        {
          std::vector<size_t> factorIds = (*it)->getFactorIndices();

          col = sslib.findCompressedIndex(sslib.getMoment("DD", popIds, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, (2. + popIdPower / 2.)));
        }

        if(popIdPower > 0)
        {
          std::vector<size_t> factorIds = (*it)->getFactorIndices();
          sslib.dropFactorIds(factorIds, id, 1);

          col = sslib.findCompressedIndex(sslib.getMoment("DD", popIds, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, -popIdPower / 2.));
        }
      }

      else if((*it)->getPrefix() == "Hl")
      {
        if(popIdCount == 2)
        {
          if((*it)->getFactorPower() < sslib.getFactorOrder())
          {
            std::vector<size_t> factorIds = (*it)->getFactorIndices();
            factorIds.push_back(id);

            col = sslib.findCompressedIndex(sslib.getMoment("Hl", popIds, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, (1. + popIdPower / 2.)));
          }

          else // truncate
          {
            std::vector<size_t> factorIds = (*it)->getFactorIndices();

            col = sslib.findCompressedIndex(sslib.getMoment("Hl", popIds, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, (1. + popIdPower / 2.)));
          }

          if(popIdPower > 0)
          {
            std::vector<size_t> factorIds = (*it)->getFactorIndices();
            sslib.dropFactorIds(factorIds, id, 1);

            col = sslib.findCompressedIndex(sslib.getMoment("Hl", popIds, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, -popIdPower / 2.));
          }
        }

        else if(popIdCount == 1)
        {
          coeffs.emplace_back(Eigen::Triplet<double>(row, row, -1. / 2.));

          if((*it)->getFactorPower() < sslib.getFactorOrder())
          {
            std::vector<size_t> factorIds = (*it)->getFactorIndices();
            factorIds.push_back(id);

            col = sslib.findCompressedIndex(sslib.getMoment("Hl", popIds, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, (1. + popIdPower) / 2.));
          }

          else // truncate
          {
            std::vector<size_t> factorIds = (*it)->getFactorIndices();

            col = sslib.findCompressedIndex(sslib.getMoment("Hl", popIds, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, (1. + popIdPower) / 2.));
          }

          if(popIdPower > 0)
          {
            std::vector<size_t> factorIds = (*it)->getFactorIndices();
            sslib.dropFactorIds(factorIds, id, 1);

            col = sslib.findCompressedIndex(sslib.getMoment("Hl", popIds, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, -popIdPower / 2.));
          }
        }

        else if(popIdCount == 0 && popIdPower > 0)
        {
          if((*it)->getFactorPower() < sslib.getFactorOrder())
          {
            std::vector<size_t> factorIds = (*it)->getFactorIndices();
            factorIds.push_back(id);

            col = sslib.findCompressedIndex(sslib.getMoment("Hl", popIds, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, popIdPower / 2.));
          }

          else // truncate
          {
            std::vector<size_t> factorIds = (*it)->getFactorIndices();

            col = sslib.findCompressedIndex(sslib.getMoment("Hl", popIds, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, popIdPower / 2.));
          }

          std::vector<size_t> factorIds = (*it)->getFactorIndices();
          sslib.dropFactorIds(factorIds, id, 1);

          col = sslib.findCompressedIndex(sslib.getMoment("Hl", popIds, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, -popIdPower / 2.));
        }
      }

      else if((*it)->getPrefix() == "Hr")
      {
        if(popIdCount == 2)
        {
          std::vector<size_t> factorIds = (*it)->getFactorIndices();

          col = sslib.findCompressedIndex(sslib.getMoment("Dr", popIds, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));
        }

        else if(popIdCount == 1)
        {
          std::vector<size_t> factorIds = (*it)->getFactorIndices();

          col = sslib.findCompressedIndex(sslib.getMoment("Dr", popIds, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1. / 2.));

          #ifdef NAKED_D
          col = sslib.findCompressedIndex(sslib.getMoment("D", { id }, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1. / 2.));
          #endif
        }
      }

      else if((*it)->getPrefix() == "pi2")
      {
        auto tmpPi2 = std::dynamic_pointer_cast<Pi2Moment>(*it);
        assert(tmpPi2 != nullptr);

        size_t countLeft = tmpPi2->getLeftHetStat()->countInstances(id);
        size_t countRight = tmpPi2->getRightHetStat()->countInstances(id);

        if(countLeft == 1) // has self-contribution, either -s/2 or +s/2
          coeffs.emplace_back(Eigen::Triplet<double>(row, row, std::pow(-1, popIds[0] != id) / 2.));

        // pi2 contributions
        if((*it)->getFactorPower() < sslib.getFactorOrder())
        {
          std::vector<size_t> factorIds = (*it)->getFactorIndices();
          factorIds.push_back(id);

          col = sslib.findCompressedIndex(sslib.getMoment("pi2", popIds, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, (countLeft + popIdPower) / 2.));
        }

        else
        {
          std::vector<size_t> factorIds = (*it)->getFactorIndices();

          col = sslib.findCompressedIndex(sslib.getMoment("pi2", popIds, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, (countLeft + popIdPower) / 2.));
        }

        if(popIdPower > 0)
        {
          std::vector<size_t> factorIds = (*it)->getFactorIndices();
          sslib.dropFactorIds(factorIds, id, 1);

          col = sslib.findCompressedIndex(sslib.getMoment("pi2", popIds, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, -popIdPower / 2.));
        }

        // contributions from other moments (D and Dr)
        // case: pi2_id_id_id_id
        if((countLeft + countRight) == 4)
        {
          // Dr contributions
          std::vector<size_t> factorIds = (*it)->getFactorIndices();

          col = sslib.findCompressedIndex(sslib.getMoment("Dr", { popIds[0], popIds[1] }, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1. / 4.));

          if((*it)->getFactorPower() < sslib.getFactorOrder())
          {
            factorIds.push_back(id);
            factorIds.push_back(id);

            col = sslib.findCompressedIndex(sslib.getMoment("Dr", { popIds[2], popIds[3] }, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1. / 4.));
          }

          else
          {
            factorIds.push_back(id);

            col = sslib.findCompressedIndex(sslib.getMoment("Dr", { popIds[2], popIds[3] }, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1. / 4.));
          }
        }

        // case: pi2_id_id_id_*
        else if(countLeft == 2 && countRight == 1)
        {
          std::vector<size_t> factorIds = (*it)->getFactorIndices();

          col = sslib.findCompressedIndex(sslib.getMoment("Dr", { popIds[2], popIds[3] }, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1. / 8.));

          if((*it)->getFactorPower() < sslib.getFactorOrder())
          {
            factorIds.push_back(id);
            factorIds.push_back(id);

            col = sslib.findCompressedIndex(sslib.getMoment("Dr", { popIds[2], popIds[3] }, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1. / 8.));
          }

          else // NOTE truncation
          {
            factorIds.push_back(id);

            col = sslib.findCompressedIndex(sslib.getMoment("Dr", { popIds[2], popIds[3] }, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1. / 8.));
          }

          /* cancel out due to right-locus permutation
          #ifdef NAKED_D
          factorIds = (*it)->getFactorIndices(); // reset

          col = sslib.findCompressedIndex(sslib.getMoment("D", { id }, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1. / 4.));

          if((*it)->getFactorPower() < sslib.getFactorOrder() - 1)
          {
            factorIds.push_back(id);
            factorIds.push_back(id);

            col = sslib.findCompressedIndex(sslib.getMoment("D", { id }, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1. / 4.));
          }

          else if((*it)->getFactorPower() < sslib.getFactorOrder())
          {
            factorIds.push_back(id);

            col = sslib.findCompressedIndex(sslib.getMoment("D", { id }, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1. / 4.));
          }

          else // NOTE heavy truncation makes us collect twice from D_id_{factorIds}
          {
            col = sslib.findCompressedIndex(sslib.getMoment("D", { id }, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1. / 4.));
          }
          #endif
          */
        }

        // case: pi2_id_*_id_id
        else if(popIds[0] == id && countRight == 2)
        {
          // Dr contributions
          std::vector<size_t> factorIds = (*it)->getFactorIndices();

          col = sslib.findCompressedIndex(sslib.getMoment("Dr", { popIds[0], popIds[1] }, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1. / 4.));

          if((*it)->getFactorPower() < sslib.getFactorOrder())
          {
            factorIds.push_back(id);
            factorIds.push_back(id);

            col = sslib.findCompressedIndex(sslib.getMoment("Dr", { popIds[2], popIds[3] }, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1. / 4.));
          }

          else
          {
            factorIds.push_back(id);

            col = sslib.findCompressedIndex(sslib.getMoment("Dr", { popIds[2], popIds[3] }, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1. / 4.));
          }
        }

        // case: pi2_id_*_id_*
        else if(popIds[0] == id && popIds[1] != id && popIds[2] == id && popIds[3] != id)
        {
          std::vector<size_t> factorIds = (*it)->getFactorIndices();

          col = sslib.findCompressedIndex(sslib.getMoment("Dr", { id, popIds[1] }, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1. / 8.));

          factorIds.push_back(id);
          col = sslib.findCompressedIndex(sslib.getMoment("Dr", { id, popIds[1] }, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1. / 8.));

          factorIds.push_back(popIds[1]);
          col = sslib.findCompressedIndex(sslib.getMoment("Dr", { id, popIds[1] }, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1. / 8.));

          sslib.dropFactorIds(factorIds, id, 1);
          col = sslib.findCompressedIndex(sslib.getMoment("Dr", { id, popIds[1] }, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1. / 8.));

          #ifdef NAKED_D
          factorIds = (*it)->getFactorIndices(); // reset

          col = sslib.findCompressedIndex(sslib.getMoment("D", { id }, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1. / 8.));

          factorIds.push_back(id);
          col = sslib.findCompressedIndex(sslib.getMoment("D", { id }, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1. / 8.));

          factorIds.push_back(popIds[1]);
          col = sslib.findCompressedIndex(sslib.getMoment("D", { id }, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1. / 8.));

          sslib.dropFactorIds(factorIds, id, 1);
          col = sslib.findCompressedIndex(sslib.getMoment("D", { id }, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1. / 8.));
          #endif
        }

        // case: pi2_*_id_id_id
        else if(popIds[0] != id && popIds[1] == id && countRight == 2)
        {
          // Dr contributions
          std::vector<size_t> factorIds = (*it)->getFactorIndices();

          col = sslib.findCompressedIndex(sslib.getMoment("Dr", { id, id }, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1. / 4.));

          factorIds.push_back(id);
          col = sslib.findCompressedIndex(sslib.getMoment("Dr", { id, id }, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1. / 4.));

          factorIds.push_back(popIds[0]);
          col = sslib.findCompressedIndex(sslib.getMoment("Dr", { id, id }, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1. / 4.));

          sslib.dropFactorIds(factorIds, id, 1);
          col = sslib.findCompressedIndex(sslib.getMoment("Dr", { id, id }, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1. / 4.));
        }

        // case: pi2_*_id_id_*
        else if(popIds[0] != id && popIds[1] == id && countRight == 1)
        {
          std::vector<size_t> factorIds = (*it)->getFactorIndices();

          col = sslib.findCompressedIndex(sslib.getMoment("Dr", { id, popIds[0] }, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1. / 8.));

          factorIds.push_back(id);
          col = sslib.findCompressedIndex(sslib.getMoment("Dr", { id, popIds[0] }, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1. / 8.));

          factorIds.push_back(popIds[1]);
          col = sslib.findCompressedIndex(sslib.getMoment("Dr", { id, popIds[0] }, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1. / 8.));

          sslib.dropFactorIds(factorIds, id, 1);
          col = sslib.findCompressedIndex(sslib.getMoment("Dr", { id, popIds[0] }, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1. / 8.));

          #ifdef NAKED_D
          factorIds = (*it)->getFactorIndices(); // reset

          col = sslib.findCompressedIndex(sslib.getMoment("D", { id }, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1. / 8.));

          factorIds.push_back(id);
          col = sslib.findCompressedIndex(sslib.getMoment("D", { id }, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1. / 8.));

          factorIds.push_back(popIds[1]);
          col = sslib.findCompressedIndex(sslib.getMoment("D", { id }, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1. / 8.));

          sslib.dropFactorIds(factorIds, id, 1);
          col = sslib.findCompressedIndex(sslib.getMoment("D", { id }, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1. / 8.));
          #endif
        }

        // case: pi2_*_*_id_id
        else if(countLeft == 0 && countRight == 2)
        {
          std::vector<size_t> factorIds = (*it)->getFactorIndices();

          col = sslib.findCompressedIndex(sslib.getMoment("Dr", { id, id }, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1. / 4.));

          factorIds.push_back(popIds[0]);
          factorIds.push_back(popIds[1]);
          col = sslib.findCompressedIndex(sslib.getMoment("Dr", { id, id }, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1. / 4.));
        }

        // case: pi2_*_*_id_*
        else if(countLeft == 0 && countRight == 1)
        {
          std::vector<size_t> factorIds = (*it)->getFactorIndices();

          col = sslib.findCompressedIndex(sslib.getMoment("Dr", { id, popIds[0] }, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1. / 8.));

          factorIds.push_back(popIds[0]);
          factorIds.push_back(popIds[1]);

          col = sslib.findCompressedIndex(sslib.getMoment("Dr", { id, popIds[0] }, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1. / 8.));

          #ifdef NAKED_D
          factorIds = (*it)->getFactorIndices(); // reset

          col = sslib.findCompressedIndex(sslib.getMoment("D", { id }, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1. / 8.));

          factorIds.push_back(popIds[0]);
          factorIds.push_back(popIds[1]);

          col = sslib.findCompressedIndex(sslib.getMoment("D", { id }, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1. / 8.));
          #endif
        }
      }

      else if((*it)->getPrefix() != "I")
        throw bpp::Exception("Selection::mis-specified Moment prefix: " + (*it)->getPrefix());
    }

    Eigen::SparseMatrix<double> mat(sizeOfBasis, sizeOfBasis);
    mat.setFromTriplets(std::begin(coeffs), std::end(coeffs));
    mat.makeCompressed();
    mat *= getParameterValue("s_" + bpp::TextTools::toString(id));
    matrices_.emplace_back(mat);
  }

  setIdentity_(sizeOfBasis);
  assembleTransitionMatrix_();
}

void Selection::updateMatrices_()
{
  for(size_t i = 0; i < matrices_.size(); ++i)
  {
    size_t id = popIndices_[i];
    std::string paramName = "s_" + bpp::TextTools::toString(id);

    double prevVal = prevParams_.getParameterValue(paramName);
    double newVal = getParameterValue(paramName);

    if(newVal != prevVal)
      matrices_[i] *= (newVal / prevVal);
  }

  assembleTransitionMatrix_();
  prevParams_.matchParametersValues(getParameters());
}
