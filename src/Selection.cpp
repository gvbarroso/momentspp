/*
 * Authors: Gustavo V. Barroso
 * Created: 22/08/2022
 * Last modified: 29/08/2023
 *
 */


#include "Selection.hpp"


// uses zero-order moment-closure approximation
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

      int power = (*it)->getPopFactorPower(id); // count of (1-2p_x) factors on focal moment

      std::vector<size_t> popIds(0);
      std::vector<size_t> factorIds = (*it)->getFactorIndices();

      if((*it)->getPrefix() == "DD")
      {
        if(power < sslib.getFactorOrder())
        {
          popIds = { (*it)->getPopIndices()[0], (*it)->getPopIndices()[1] };
          factorIds.push_back(id);

          col = sslib.findCompressedIndex(sslib.getMoment("DD", popIds, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, (2. + power/2.)));
        }

        else
        {
          popIds = { (*it)->getPopIndices()[0], (*it)->getPopIndices()[1] };

          col = sslib.findCompressedIndex(sslib.getMoment("DD", popIds, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, (2. + power/2.)));
        }

        if(power > 0)
        {
          popIds = { (*it)->getPopIndices()[0], (*it)->getPopIndices()[1] };
          factorIds = (*it)->getFactorIndices();
          sslib.dropFactorIds(factorIds, id, 1);

          col = sslib.findCompressedIndex(sslib.getMoment("DD", popIds, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, -power/2.));
        }
      }

      else if((*it)->getPrefix() == "Dr")
      {
        if(power < sslib.getFactorOrder())
        {
          popIds = { (*it)->getPopIndices()[0], (*it)->getPopIndices()[1] };
          factorIds.push_back(id);

          col = sslib.findCompressedIndex(sslib.getMoment("Dr", popIds, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, (1. + power/2.)));
        }

        else
        {
          popIds = { (*it)->getPopIndices()[0], (*it)->getPopIndices()[1] };

          col = sslib.findCompressedIndex(sslib.getMoment("Dr", popIds, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, (1. + power/2.)));
        }

        if(power > 0)
        {
          popIds = { (*it)->getPopIndices()[0], (*it)->getPopIndices()[1] };
          factorIds = (*it)->getFactorIndices();
          sslib.dropFactorIds(factorIds, id, 1);

          col = sslib.findCompressedIndex(sslib.getMoment("Dr", popIds, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, -power/2.));
        }

        if(power < sslib.getFactorOrder())
        {
          popIds = { (*it)->getPopIndices()[0], (*it)->getPopIndices()[1] };
          factorIds = (*it)->getFactorIndices();

          col = sslib.findCompressedIndex(sslib.getMoment("DD", popIds, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, -2.));
        }

        else
        {
          popIds = { (*it)->getPopIndices()[0], (*it)->getPopIndices()[1] };
          factorIds = (*it)->getFactorIndices();
          sslib.dropFactorIds(factorIds, id, 1);

          col = sslib.findCompressedIndex(sslib.getMoment("DD", popIds, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, -2.));
        }
      }

      else if((*it)->getPrefix() == "Hl")
      {
        if(power < sslib.getFactorOrder())
        {
          popIds = { (*it)->getPopIndices()[0], (*it)->getPopIndices()[1] };
          factorIds.push_back(id);

          col = sslib.findCompressedIndex(sslib.getMoment("Hl", popIds, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, (1. + power/2.)));
        }

        else
        {
          popIds = { (*it)->getPopIndices()[0], (*it)->getPopIndices()[1] };

          col = sslib.findCompressedIndex(sslib.getMoment("Hl", popIds, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, (1. + power/2.)));
        }

        if(power > 0)
        {
          popIds = { (*it)->getPopIndices()[0], (*it)->getPopIndices()[1] };
          factorIds = (*it)->getFactorIndices();
          sslib.dropFactorIds(factorIds, id, 1);

          col = sslib.findCompressedIndex(sslib.getMoment("Hl", popIds, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, -power/2.));
        }
      }

      else if((*it)->getPrefix() == "Hr")
      {
        popIds = { (*it)->getPopIndices()[0], (*it)->getPopIndices()[1] };

        col = sslib.findCompressedIndex(sslib.getMoment("Dr", popIds, { 0 }));
        coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));
      }

      else if((*it)->getPrefix() == "pi2")
      {
        if(power < sslib.getFactorOrder())
        {
          popIds = { (*it)->getPopIndices()[0], (*it)->getPopIndices()[1], (*it)->getPopIndices()[2], (*it)->getPopIndices()[3] };
          factorIds.push_back(id);

          col = sslib.findCompressedIndex(sslib.getMoment("pi2", popIds, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, (1. + power/2.)));
        }

        else
        {
          popIds = { (*it)->getPopIndices()[0], (*it)->getPopIndices()[1], (*it)->getPopIndices()[2], (*it)->getPopIndices()[3] };

          col = sslib.findCompressedIndex(sslib.getMoment("pi2", popIds, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, (1. + power/2.)));
        }

        if(power > 0)
        {
          popIds = { (*it)->getPopIndices()[0], (*it)->getPopIndices()[1], (*it)->getPopIndices()[2], (*it)->getPopIndices()[3] };
          factorIds = (*it)->getFactorIndices();
          sslib.dropFactorIds(factorIds, id, 1);

          col = sslib.findCompressedIndex(sslib.getMoment("pi2", popIds, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, -power/2.));
        }

        // WARNING on contributions from Dr collecting only from pop ids 0 and 1 (left), maybe makes sense?
        popIds = { (*it)->getPopIndices()[0], (*it)->getPopIndices()[1] };
        factorIds = (*it)->getFactorIndices();

        col = sslib.findCompressedIndex(sslib.getMoment("Dr", popIds, factorIds));
        coeffs.emplace_back(Eigen::Triplet<double>(row, col, 0.25));

        if(power < sslib.getFactorOrder())
        {
          popIds = { (*it)->getPopIndices()[0], (*it)->getPopIndices()[1] };
          factorIds.push_back(id);
          factorIds.push_back(id);

          col = sslib.findCompressedIndex(sslib.getMoment("Dr", popIds, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, -0.25));
        }

        else
        {
          popIds = { (*it)->getPopIndices()[0], (*it)->getPopIndices()[1] };
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
