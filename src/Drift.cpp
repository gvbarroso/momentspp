/*
 * Authors: Gustavo V. Barroso
 * Created: 09/08/2022
 * Last modified: 21/04/2023
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

    for(auto it = std::begin(sslib.getBasis()); it != std::end(sslib.getBasis()); ++it)
    {
      int row = it - std::begin(sslib.getBasis());
      int col = -1;
      size_t popIdCount = (*it)->countInstances(id); // count of i (focal pop ID) in moment's name

      if((*it)->getPrefix() == "DD")
      {
        if(popIdCount == 2)
        {
          coeffs.emplace_back(Eigen::Triplet<double>(row, row, -3.));

          col = sslib.findCompressedIndex(sslib.findDzIndex(id, id, id));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));

          col = sslib.findCompressedIndex(sslib.findPi2Index(id, id, id, id));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));
        }

        else if(popIdCount == 1)
          coeffs.emplace_back(Eigen::Triplet<double>(row, row, -1.));
      }

      else if((*it)->getPrefix() == "Dz")
      {
        if((*it)->getPopIndices()[0] == id) // if D_i_z**
        {
          if(popIdCount == 3) // D_i_z_ii
          {
            coeffs.emplace_back(Eigen::Triplet<double>(row, row, -5.));

            col = sslib.findCompressedIndex(sslib.findDdIndex(id, id));
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
            col = sslib.findCompressedIndex(sslib.findDdIndex(id, (*it)->getPopIndices()[0]));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, 2.));

            col = sslib.findCompressedIndex(sslib.findDdIndex((*it)->getPopIndices()[0], id));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, 2.));
          }
        }
      }

      else if((*it)->getPrefix() == "pi2")
      {
        auto tmpPi2 = std::dynamic_pointer_cast<Pi2Moment>(*it);
        assert(tmpPi2 != nullptr);

        size_t countLeft = tmpPi2->getLeftHetStat()->countInstances(id);
        size_t countRight = tmpPi2->getRightHetStat()->countInstances(id);

        if((countLeft + countRight) == 4)
        {
          coeffs.emplace_back(Eigen::Triplet<double>(row, row, -2.));

          col = sslib.findCompressedIndex(sslib.findDzIndex(id, id, id));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));
        }

        else if((countLeft == 2) || (countRight == 2))
        {
          coeffs.emplace_back(Eigen::Triplet<double>(row, row, -1.));

          if((countLeft + countRight) == 3)
          {
            for(size_t j = 0; j < numPops; ++j)
            {
              size_t jd = popIndices_[j];

              if(id != jd && (*it)->hasPopIndex(jd))
              {
                col = sslib.findCompressedIndex(sslib.findDzIndex(id, id, jd));
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, 0.25));

                col = sslib.findCompressedIndex(sslib.findDzIndex(id, jd, id));
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, 0.25));
              }
            }
          }
        }

        else if(countLeft == 1 && countRight == 1)
        {
          for(size_t j = 0; j < numPops; ++j)
          {
            size_t jd = popIndices_[j];

            if(id != jd && (*it)->hasPopIndex(jd))
            {
              std::vector<size_t> diff = (*it)->fetchDiffPopIds(id);
              assert(diff.size() == 2);
              double f = 1. + (diff[0] == diff[1]);

              col = sslib.findCompressedIndex(sslib.findDzIndex(id, diff[0], diff[1]));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, f * 0.0625));

              col = sslib.findCompressedIndex(sslib.findDzIndex(id, diff[1], diff[0]));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, f * 0.0625));
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
