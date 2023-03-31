/*
 * Authors: Gustavo V. Barroso
 * Created: 09/08/2022
 * Last modified: 21/03/2023
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
    std::vector<Eigen::Triplet<double>> coeffs(0);
    coeffs.reserve(sizeOfBasis);

    for(auto it = std::begin(sslib.getBasis()); it != std::end(sslib.getBasis()); ++it)
    {
      int row = it - std::begin(sslib.getBasis());
      int col = -1;
      size_t popIdCount = (*it)->countInstances(i); // count of i (focal pop ID) in moment's name

      if((*it)->getPrefix() == "DD")
      {
        if(popIdCount == 2)
        {
          coeffs.emplace_back(Eigen::Triplet<double>(row, row, -3.));

          col = sslib.findCompressedIndex(sslib.findDzIndex(i, i, i));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));

          col = sslib.findCompressedIndex(sslib.findPi2Index(i, i, i, i));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));
        }

        else if(popIdCount == 1)
          coeffs.emplace_back(Eigen::Triplet<double>(row, row, -1.));
      }

      else if((*it)->getPrefix() == "Dz")
      {
        if((*it)->getPopIndices()[0] == i) // if D_i_z**
        {
          if(popIdCount == 3) // D_i_z_ii
          {
            coeffs.emplace_back(Eigen::Triplet<double>(row, row, -5.));

            col = sslib.findCompressedIndex(sslib.findDdIndex(i, i));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, 4.));
          }

          else if(popIdCount == 2) // D_i_z_i* or D_i_z_*i
            coeffs.emplace_back(Eigen::Triplet<double>(row, row, -3.));

          else if(popIdCount == 1) // if D_i_z_xx
            coeffs.emplace_back(Eigen::Triplet<double>(row, row, -1.));
        }

        else // if D_x_z**
        {
          if(popIdCount == 2) // D_x_z_ii
          {
            col = sslib.findCompressedIndex(sslib.findDdIndex(i, (*it)->getPopIndices()[0]));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, 2.));

            col = sslib.findCompressedIndex(sslib.findDdIndex((*it)->getPopIndices()[0], i));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, 2.));
          }
        }
      }

      else if((*it)->getPrefix() == "pi2")
      {
        auto tmpPi2 = std::dynamic_pointer_cast<Pi2Moment>(*it);
        assert(tmpPi2 != nullptr);

        size_t countLeft = tmpPi2->getLeftHetStat()->countInstances(i);
        size_t countRight = tmpPi2->getRightHetStat()->countInstances(i);

        if((countLeft + countRight) == 4)
        {
          coeffs.emplace_back(Eigen::Triplet<double>(row, row, -2.));

          col = sslib.findCompressedIndex(sslib.findDzIndex(i, i, i));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));
        }

        else if((countLeft == 2) || (countRight == 2))
        {
          coeffs.emplace_back(Eigen::Triplet<double>(row, row, -1.));

          if((countLeft + countRight) == 3)
          {
            for(size_t j = 0; j < numPops; ++j)
            {
              if(i != j)
              {
                col = sslib.findCompressedIndex(sslib.findDzIndex(i, i, j));
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, 0.25));

                col = sslib.findCompressedIndex(sslib.findDzIndex(i, j, i));
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, 0.25));
              }
            }
          }
        }

        else if(countLeft == 1 && countRight == 1)
        {
          for(size_t j = 0; j < numPops; ++j)
          {
            if(i != j)
            {
              col = sslib.findCompressedIndex(sslib.findDzIndex(i, j, j));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, 0.25));
            }
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
    mat *= (getParameterValue("1/2N_" + bpp::TextTools::toString(i)));
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
    paramName = "1/2N_" + bpp::TextTools::toString(i);

    double prevVal = prevParams_.getParameterValue(paramName);
    double newVal = getParameterValue(paramName);

    if(newVal != prevVal)
      matrices_[i] *= (newVal / prevVal);
  }

  assembleTransitionMatrix_();
  prevParams_.matchParametersValues(getParameters());
}
