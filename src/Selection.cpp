/*
 * Authors: Gustavo V. Barroso
 * Created: 22/08/2022
 * Last modified: 07/09/2023
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
      std::vector<size_t> factorIds = (*it)->getFactorIndices();

      if((*it)->getPrefix() == "DD")
      {
        if((*it)->getFactorPower() < sslib.getFactorOrder())
        {
          factorIds.push_back(id);

          col = sslib.findCompressedIndex(sslib.getMoment("DD", popIds, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, (popIdCount + popIdPower/2.)));
        }

        else // truncate
        {
          col = sslib.findCompressedIndex(sslib.getMoment("DD", popIds, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, (2. + popIdPower/2.)));
        }

        if(popIdPower > 0)
        {
          factorIds = (*it)->getFactorIndices(); // reset
          sslib.dropFactorIds(factorIds, id, 1);

          col = sslib.findCompressedIndex(sslib.getMoment("DD", popIds, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, -popIdPower/2.));
        }
      }

      else if((*it)->getPrefix() == "Dr")
      {
        if(popIdCount == 2)
        {
          // Dr contributions
          if((*it)->getFactorPower() < sslib.getFactorOrder())
          {
            factorIds.push_back(id);

            col = sslib.findCompressedIndex(sslib.getMoment("Dr", popIds, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, (1. + popIdPower/2.)));
          }

          else // truncate
          {
            col = sslib.findCompressedIndex(sslib.getMoment("Dr", popIds, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, (1. + popIdPower/2.)));
          }

          if(popIdPower > 0)
          {
            factorIds = (*it)->getFactorIndices(); // reset
            sslib.dropFactorIds(factorIds, id, 1);

            col = sslib.findCompressedIndex(sslib.getMoment("Dr", popIds, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, -popIdPower/2.));
          }

          // DD contributions
          if((*it)->getFactorPower() < sslib.getFactorOrder())
          {
            factorIds = (*it)->getFactorIndices(); // reset

            col = sslib.findCompressedIndex(sslib.getMoment("DD", popIds, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, -2.));
          }

          else if(popIdPower > 0) // NOTE used to be just 'else'
          {
            factorIds = (*it)->getFactorIndices(); // reset
            sslib.dropFactorIds(factorIds, id, 1);

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
              factorIds.push_back(id);

              col = sslib.findCompressedIndex(sslib.getMoment("Dr", popIds, factorIds));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, (1. + popIdPower/2.)));
            }

            else // truncate
            {
              col = sslib.findCompressedIndex(sslib.getMoment("Dr", popIds, factorIds));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, (1. + popIdPower/2.)));
            }

            if(popIdPower > 0)
            {
              factorIds = (*it)->getFactorIndices(); // reset
              sslib.dropFactorIds(factorIds, id, 1);

              col = sslib.findCompressedIndex(sslib.getMoment("Dr", popIds, factorIds));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, -popIdPower/2.));
            }
          }

          else // if(popIds[1] == id)
          {
            if((*it)->getFactorPower() < sslib.getFactorOrder())
            {
              factorIds.push_back(id);

              col = sslib.findCompressedIndex(sslib.getMoment("Dr", popIds, factorIds));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, (1. + popIdPower/2.)));
            }

            else // truncate
            {
              col = sslib.findCompressedIndex(sslib.getMoment("Dr", popIds, factorIds));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, (1. + popIdPower/2.)));
            }

            if(popIdPower > 0)
            {
              factorIds = (*it)->getFactorIndices(); // reset
              sslib.dropFactorIds(factorIds, id, 1);

              col = sslib.findCompressedIndex(sslib.getMoment("Dr", popIds, factorIds));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, -popIdPower/2.));
            }

            // DD contributions
            if((*it)->getFactorPower() < sslib.getFactorOrder())
            {
              factorIds = (*it)->getFactorIndices(); // reset

              col = sslib.findCompressedIndex(sslib.getMoment("DD", popIds, factorIds));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, -2.));
            }

            else if(popIdPower > 0) // NOTE used to be just 'else'
            {
              factorIds = (*it)->getFactorIndices(); // reset
              sslib.dropFactorIds(factorIds, id, 1);

              col = sslib.findCompressedIndex(sslib.getMoment("DD", popIds, factorIds));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, -2.));
            }
          }
        }
      }

      else if((*it)->getPrefix() == "Hl")
      {
        if(popIdCount == 2)
        {
          if((*it)->getFactorPower() < sslib.getFactorOrder())
          {
            factorIds.push_back(id);

            col = sslib.findCompressedIndex(sslib.getMoment("Hl", popIds, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, (1. + popIdPower/2.)));
          }

          else // truncate
          {
            col = sslib.findCompressedIndex(sslib.getMoment("Hl", popIds, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, (1. + popIdPower/2.)));
          }

          if(popIdPower > 0)
          {
            factorIds = (*it)->getFactorIndices(); // reset
            sslib.dropFactorIds(factorIds, id, 1);

            col = sslib.findCompressedIndex(sslib.getMoment("Hl", popIds, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, -popIdPower/2.));
          }
        }

        else
        {
          coeffs.emplace_back(Eigen::Triplet<double>(row, row, -1./2.));

          if((*it)->getFactorPower() < sslib.getFactorOrder())
          {
            factorIds.push_back(id);

            col = sslib.findCompressedIndex(sslib.getMoment("Hl", popIds, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, (1. + popIdPower)/2.));
          }

          else // truncate
          {
            col = sslib.findCompressedIndex(sslib.getMoment("Hl", popIds, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, (1. + popIdPower)/2.));
          }

          if(popIdPower > 0)
          {
            factorIds = (*it)->getFactorIndices(); // reset
            sslib.dropFactorIds(factorIds, id, 1);

            col = sslib.findCompressedIndex(sslib.getMoment("Hl", popIds, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, -popIdPower/2.));
          }
        }
      }

      else if((*it)->getPrefix() == "Hr")
      {
        if(popIdCount == 2)
        {
          col = sslib.findCompressedIndex(sslib.getMoment("Dr", popIds, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));
        }

        else if(popIdCount == 1)
        {
          col = sslib.findCompressedIndex(sslib.getMoment("Dr", popIds, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./2.));

          col = sslib.findCompressedIndex(sslib.getMoment("D", { id }, factorIds)); // TODO include D^1 stats in Basis
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./2.));

          if(popIdPower > 0)
          {
            if((*it)->getFactorPower() < sslib.getFactorOrder())
            {
              factorIds.push_back(id);

              col = sslib.findCompressedIndex(sslib.getMoment("Hr", popIds, factorIds));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, (popIdPower/2.)));
            }

            else // truncate
            {
              col = sslib.findCompressedIndex(sslib.getMoment("Hr", popIds, factorIds));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, (popIdPower/2.)));
            }

            factorIds = (*it)->getFactorIndices(); // reset
            sslib.dropFactorIds(factorIds, id, 1);

            col = sslib.findCompressedIndex(sslib.getMoment("Hr", popIds, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, -popIdPower/2.));
          }
        }
      }

      else if((*it)->getPrefix() == "pi2")
      {
        if((*it)->getFactorPower() < sslib.getFactorOrder())
        {
          factorIds.push_back(id);

          col = sslib.findCompressedIndex(sslib.getMoment("pi2", popIds, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, (1. + popIdPower/2.)));
        }

        else
        {
          col = sslib.findCompressedIndex(sslib.getMoment("pi2", popIds, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, (1. + popIdPower/2.)));
        }

        if(popIdPower > 0)
        {
          factorIds = (*it)->getFactorIndices();
          sslib.dropFactorIds(factorIds, id, 1);

          col = sslib.findCompressedIndex(sslib.getMoment("pi2", popIds, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, -popIdPower/2.));
        }

        // NOTE on contributions from Dr collecting only from pop ids 0 and 1 (left), maybe makes sense?
        factorIds = (*it)->getFactorIndices();

        col = sslib.findCompressedIndex(sslib.getMoment("Dr", popIds, factorIds));
        coeffs.emplace_back(Eigen::Triplet<double>(row, col, 0.25));

        if((*it)->getFactorPower() < sslib.getFactorOrder())
        {
          factorIds.push_back(id);
          factorIds.push_back(id);

          col = sslib.findCompressedIndex(sslib.getMoment("Dr", popIds, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, -0.25));
        }

        else
        {
          factorIds.push_back(id);

          col = sslib.findCompressedIndex(sslib.getMoment("Dr", popIds, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, -0.25));
        }
      }

      else if((*it)->getPrefix() != "I" && (*it)->getPrefix() != "DD" && (*it)->getPrefix() != "Dr")
        throw bpp::Exception("Mutation::mis-specified Moment prefix: " + (*it)->getPrefix());
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
  std::string paramName = "s";

  double prevVal = prevParams_.getParameterValue(paramName);
  double newVal = getParameterValue(paramName);

  matrices_.front() *= (newVal / prevVal);
  prevParams_.matchParametersValues(getParameters());
}
