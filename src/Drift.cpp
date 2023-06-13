/*
 * Authors: Gustavo V. Barroso
 * Created: 09/08/2022
 * Last modified: 13/06/2023
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

    int j = 0;
    int k = -3;

    int l = -1;
    int m = -2;

    int n = -1;
    int o = -1;

    for(auto it = std::begin(sslib.getBasis()); it != std::end(sslib.getBasis()); ++it)
    {
      int row = it - std::begin(sslib.getBasis());
      int col = -1;

      size_t popIdCount = (*it)->countInstances(id); // count of i (focal pop ID) in moment's name
      size_t x = (*it)->getFactorPower(); // count of (1-2p) factors on focal moment

      if((*it)->getPrefix() == "DD")
      {
        if(popIdCount == 2)
        {
          j += k;
          --k;

          coeffs.emplace_back(Eigen::Triplet<double>(row, row, j));

          if(x >= 2)
          {
            double y = (x * (x - 1)) / 2.;
            col = sslib.findCompressedIndex(sslib.findDdIndex(id, id, x - 2));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, y));
          }

          col = sslib.findCompressedIndex(sslib.findDrIndex(id, id, x + 1));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));

          col = sslib.findCompressedIndex(sslib.findPi2Index(id, id, id, id, x));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));
        }
      }

      else if((*it)->getPrefix() == "Dr")
      {
        if((*it)->getPopIndices()[0] == id) // D_i_r*
        {
          if(popIdCount == 2) // D_i_r_i
          {
            l += m;
            --m;

            coeffs.emplace_back(Eigen::Triplet<double>(row, row, l));

            if(x >= 1)
            {
              col = sslib.findCompressedIndex(sslib.findDdIndex(id, id, x - 1));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, 2. * x));

              if(x >= 2)
              {
                double y = (x * (x - 1)) / 2.;
                col = sslib.findCompressedIndex(sslib.findDrIndex(id, id, x - 2));
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

        if(x >= 1)
        {
          col = sslib.findCompressedIndex(sslib.findDrIndex(id, id, 1));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));

          if(x >= 2)
          {
            double y = (x * (x - 1)) / 2.;
            col = sslib.findCompressedIndex(sslib.findPi2Index(id, id, id, id, x - 2));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, y));
          }
        }
      }

      else if((*it)->getPrefix() == "H")
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
  std::string paramName = "";

  for(size_t i = 0; i < matrices_.size(); ++i)
  {
    size_t id = popIndices_[i];
    paramName = "1/2N_" + bpp::TextTools::toString(id);

    double prevVal = prevParams_.getParameterValue(paramName);
    double newVal = getParameterValue(paramName);

    if(newVal != prevVal)
      matrices_[i] *= (newVal / prevVal);
  }

  assembleTransitionMatrix_();
  prevParams_.matchParametersValues(getParameters());
}
