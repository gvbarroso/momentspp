/*
 * Authors: Gustavo V. Barroso
 * Created: 09/08/2022
 * Last modified: 28/08/2023
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

    int l = -2;
    int m = -1;

    int n = -1;
    int o = -1;

    for(auto it = std::begin(sslib.getBasis()); it != std::end(sslib.getBasis()); ++it)
    {
      int row = it - std::begin(sslib.getBasis());
      int col = -1;

      size_t popIdCount = (*it)->countInstances(id);
      int power = (*it)->getFactorPower(); // TODO getPopFactorPower() ?

      std::vector<size_t> popIds(0);
      std::vector<size_t> factorIds(0);

      if((*it)->getPrefix() == "DD")
      {
        if(popIdCount == 2)
        {
          j += k;
          --k;

          coeffs.emplace_back(Eigen::Triplet<double>(row, row, j));

          popIds = { id, id };
          col = sslib.findCompressedIndex(sslib.findMoment("Dr", popIds, power + 1));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));

          popIds = { id, id, id, id };
          col = sslib.findCompressedIndex(sslib.findMoment("pi2", popIds, power));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));

          if(power > 1)
          {
            popIds = { id, id };
            double y = (power * (power - 1)) / 2.;
            col = sslib.findCompressedIndex(sslib.findMoment("DD", popIds, power - 2));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, y));
          }
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

            if(power > 0)
            {
              popIds = { id, id };
              col = sslib.findCompressedIndex(sslib.findMoment("DD", popIds, power - 1));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, 4. * power));

              if(power > 1)
              {
                double y = (power * (power - 1)) / 2.;
                col = sslib.findCompressedIndex(sslib.findMoment("Dr", popIds, power - 2));
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
        col = sslib.findCompressedIndex(sslib.findMoment("Dr", popIds, power + 1));
        coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1. + power/2.));

        if(power > 0)
        {
          col = sslib.findCompressedIndex(sslib.findMoment("Dr", popIds, power - 1));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, (-2. * power) / 4.));

          if(power > 1)
          {
            popIds = { id, id, id, id };
            double y = (power * (power - 1)) / 2.;
            col = sslib.findCompressedIndex(sslib.findMoment("pi2", popIds, power - 2));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, y));
          }
        }
      }

      else if((*it)->getPrefix() == "Hl")
      {
        if(popIdCount == 2)
        {
          double y = - ((power + 2) * (power + 1)) / 2.;
          coeffs.emplace_back(Eigen::Triplet<double>(row, row, y));

          if(power > 1)
          {
            popIds = { id, id };
            double z = (power * (power - 1)) / 2.;
            col = sslib.findCompressedIndex(sslib.findMoment("Hl", popIds, power - 2));
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
