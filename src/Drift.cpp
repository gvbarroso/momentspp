/*
 * Authors: Gustavo V. Barroso
 * Created: 09/08/2022
 * Last modified: 25/09/2023
 *
 */

#include <Bpp/Text/TextTools.h>

#include "Drift.hpp"

void Drift::setUpMatrices_(const SumStatsLibrary& sslib)
{
  size_t numPops = getParameters().size();
  size_t sizeOfBasis = sslib.getSizeOfBasis();
  matrices_.reserve(numPops);

  for(size_t i = 0; i < numPops; ++i)
  {
    size_t id = popIndices_[i];
    std::vector<Eigen::Triplet<double>> coeffs(0);
    coeffs.reserve(sizeOfBasis);

    // helpers for counting coefficients (assumes basis is ordered with SumStatsLibrary::compareMoments_())
    // D
    int a = -1;
    int b = -2;
    int c = 1;
    int d = 2;

    // DD
    int e = -2;
    int f = -1;

    int j = 0;
    int k = -3;

    int l = -2;
    int m = -1;

    int n = -1;
    int o = -1;

    for(auto it = std::begin(sslib.getBasis()); it != std::end(sslib.getBasis()); ++it)
    {
      int row = it - std::begin(sslib.getBasis());
      int col = -1;

      size_t popIdCount = (*it)->countInstances(id);
      int popIdPower = (*it)->getPopFactorPower(id);

      std::vector<size_t> popIds(0);

      if((*it)->getPrefix() == "D")
      {
        if(popIdCount == 1)
        {
          a += b;
          --b;

          coeffs.emplace_back(Eigen::Triplet<double>(row, row, a));

          if(popIdPower > 1)
          {
            c += d;
            ++d;

            std::vector<size_t> factorIds = (*it)->getFactorIndices();
            sslib.dropFactorIds(factorIds, id, 2);

            col = sslib.findCompressedIndex(sslib.getMoment("D", { id }, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, c));
          }
        }
      }

      else if((*it)->getPrefix() == "DD")
      {
        if(popIdCount == 2)
        {
          j += k;
          --k;

          coeffs.emplace_back(Eigen::Triplet<double>(row, row, j));

          std::vector<size_t> factorIds = (*it)->getFactorIndices();
          popIds = { id, id, id, id };
          col = sslib.findCompressedIndex(sslib.getMoment("pi2", popIds, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));

          popIds = { id, id };
          factorIds.push_back(id);
          col = sslib.findCompressedIndex(sslib.getMoment("Dr", popIds, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));
        }

        else if(popIdCount == 1)
        {
          e += f;
          --f;

          coeffs.emplace_back(Eigen::Triplet<double>(row, row, e));
        }

        if(popIdPower > 1)
        {
          std::vector<size_t> factorIds = (*it)->getFactorIndices();
          sslib.dropFactorIds(factorIds, id, 2);

          double y = (popIdPower * (popIdPower - 1)) / 2.;
          col = sslib.findCompressedIndex(sslib.getMoment("DD", popIds, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, y));
        }
      }

      else if((*it)->getPrefix() == "Dr") // TODO improve
      {
        if((*it)->getPopIndices()[0] == id) // D_id_r_-
        {
          if(popIdCount == 2) // D_id_r_id
          {
            l += m;
            --m;

            coeffs.emplace_back(Eigen::Triplet<double>(row, row, l));

            if(popIdPower > 0)
            {
              popIds = { id, id };
              std::vector<size_t> factorIds = (*it)->getFactorIndices();
              sslib.dropFactorIds(factorIds, id, 1);

              col = sslib.findCompressedIndex(sslib.getMoment("DD", popIds, factorIds));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, 4. * popIdPower));

              if(popIdPower > 1)
              {
                sslib.dropFactorIds(factorIds, id, 1);

                double y = (popIdPower * (popIdPower - 1)) / 2.;
                col = sslib.findCompressedIndex(sslib.getMoment("Dr", popIds, factorIds));
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, y));
              }
            }
          }
        }
      }

      else if((*it)->getPrefix() == "pi2")
      {
        n += o;
        --o;

        coeffs.emplace_back(Eigen::Triplet<double>(row, row, n));

        popIds = { id, id };
        std::vector<size_t> factorIds = (*it)->getFactorIndices();
        factorIds.push_back(id);

        col = sslib.findCompressedIndex(sslib.getMoment("Dr", popIds, factorIds));
        coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1. + popIdPower/2.));

        if(popIdPower > 0)
        {
          sslib.dropFactorIds(factorIds, id, 2);
          col = sslib.findCompressedIndex(sslib.getMoment("Dr", popIds, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, (-2. * popIdPower) / 4.));

          if(popIdPower > 1)
          {
            popIds = { id, id, id, id };
            sslib.dropFactorIds(factorIds, id, 1);

            double y = (popIdPower * (popIdPower - 1)) / 2.;
            col = sslib.findCompressedIndex(sslib.getMoment("pi2", popIds, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, y));
          }
        }
      }

      else if((*it)->getPrefix() == "Hl")
      {
        if(popIdCount == 2)
        {
          double y = - ((popIdPower + 2) * (popIdPower + 1)) / 2.;
          coeffs.emplace_back(Eigen::Triplet<double>(row, row, y));

          if(popIdPower > 1)
          {
            popIds = { id, id };
            std::vector<size_t> factorIds = (*it)->getFactorIndices();
            sslib.dropFactorIds(factorIds, id, 2);

            double z = (popIdPower * (popIdPower - 1)) / 2.;
            col = sslib.findCompressedIndex(sslib.getMoment("Hl", popIds, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, z));
          }
        }
      }

      else if((*it)->getPrefix() == "Hr")
      {
        if(popIdCount == 2)
          coeffs.emplace_back(Eigen::Triplet<double>(row, row, -1.));
      }

      else if((*it)->getPrefix() != "I")
        throw bpp::Exception("Drift::mis-specified Moment prefix: " + (*it)->getPrefix());
    }

    Eigen::SparseMatrix<double> mat(sizeOfBasis, sizeOfBasis);
    mat.setFromTriplets(std::begin(coeffs), std::end(coeffs));
    mat.makeCompressed();
    mat *= (getParameterValue("1/2N_" + bpp::TextTools::toString(id)));
    matrices_.emplace_back(mat);
  }

  setIdentity_(sizeOfBasis);
  assembleTransitionMatrix_();
}

void Drift::updateMatrices_()
{
  for(size_t i = 0; i < matrices_.size(); ++i)
  {
    size_t id = popIndices_[i];
    std::string paramName = "1/2N_" + bpp::TextTools::toString(id);

    double prevVal = prevParams_.getParameterValue(paramName);
    double newVal = getParameterValue(paramName);

    if(newVal != prevVal)
      matrices_[i] *= (newVal / prevVal);
  }

  assembleTransitionMatrix_();
  prevParams_.matchParametersValues(getParameters());
}
