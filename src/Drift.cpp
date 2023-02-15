/*
 * Authors: Gustavo V. Barroso
 * Created: 09/08/2022
 * Last modified: 15/02/2023
 *
 */

#include <Bpp/Text/TextTools.h>

#include "Drift.hpp"


void Drift::setUpMatrices_(const SumStatsLibrary& sslib)
{
  size_t numPops = getParameters().size();
  size_t numStats = sslib.getNumStats();

  matrices_.reserve(numPops);

  // for each focal population (and cross-population terms involving it)
  for(size_t i = 0; i < numPops; ++i)
  {
    std::vector<Eigen::Triplet<double>> coeffs(0);
    coeffs.reserve(numStats);

    // for each stat in vector Y (going by rows of matrices_[i])
    for(auto it = std::begin(sslib.getMoments()); it != std::end(sslib.getMoments()); ++it)
    {
      int row = it - std::begin(sslib.getMoments()); // row index
      int col = -1; // inits column index to out-of-bounds
      size_t popIdCount = (*it)->countInstances(i); // count of i (focal pop ID) in moment's name

      if((*it)->getPrefix() == "DD")
      {
        if(popIdCount == 2)
        {
          coeffs.emplace_back(Eigen::Triplet<double>(row, row, -3.));

          col = sslib.findDzIndex(i, i, i);
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));

          // Pi2Moment permutations:
          col = sslib.findPi2Index(i, i, i, i, "A");
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));

          col = sslib.findPi2Index(i, i, i, i, "B");
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));

          col = sslib.findPi2Index(i, i, i, i, "C");
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));

          col = sslib.findPi2Index(i, i, i, i, "D");
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

            col = sslib.findDdIndex(i, i);
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
            col = sslib.findDdIndex(i, (*it)->getPopIndices()[0]);
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, 4.));

            col = sslib.findDdIndex((*it)->getPopIndices()[0], i);
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, 4.));
          }
        }
      }

      else if((*it)->getPrefix() == "pi2")
      {
        size_t countLeft = 0; // count of i before ';' character in pi2(**;**)
        size_t countRight = 0; // count of i after ';' character in pi2(**;**)

        if((*it)->getPopIndices()[0] == i)
          ++countLeft;

        if((*it)->getPopIndices()[1] == i)
          ++countLeft;

        if((*it)->getPopIndices()[2] == i)
          ++countRight;

        if((*it)->getPopIndices()[3] == i)
          ++countRight;

        if((countLeft + countRight) == 4)
        {
          coeffs.emplace_back(Eigen::Triplet<double>(row, row, -2.));

          col = sslib.findDzIndex(i, i, i);
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./4.)); // (1 / # of Pi2(iiii)_*)
        }

        else if((countLeft == 2) || (countRight == 2))
          coeffs.emplace_back(Eigen::Triplet<double>(row, row, -1.));
      }

      else if((*it)->getPrefix() == "H")
      {
        if(popIdCount == 2)
          coeffs.emplace_back(Eigen::Triplet<double>(row, row, -1.));
      }

      else if((*it)->getPrefix() != "I")
        throw bpp::Exception("Drift::mis-specified Moment prefix: " + (*it)->getPrefix());
    }

    Eigen::SparseMatrix<double> mat(numStats, numStats);
    mat.setFromTriplets(std::begin(coeffs), std::end(coeffs));
    mat.makeCompressed();
    mat *= (getParameterValue("1/2N_" + bpp::TextTools::toString(i))); // scales because of updateMatrices_()
    matrices_.emplace_back(mat); // at the i-th position of vector, where i index the population
  }
}

void Drift::updateMatrices_()
{
  std::string paramName = "";

  for(size_t i = 0; i < matrices_.size(); ++i) // one matrix / eigensolver per population (within each Epoch)
  {
    paramName = "1/2N_" + bpp::TextTools::toString(i);

    double prevVal = prevParams_.getParameterValue(paramName); // old
    double newVal = getParameterValue(paramName); // new

    matrices_[i] *= (newVal / prevVal);
  }

  prevParams_.matchParametersValues(getParameters());
}
