/*
 * Authors: Gustavo V. Barroso
 * Created: 09/08/2022
 * Last modified: 01/09/2022
 *
 */

#include <Bpp/Text/TextTools.h>

#include "Drift.hpp"


void Drift::setUpMatrices_(const SumStatsLibrary& sslib)
{
  size_t numPops = sslib.getNumPops();

  matrices_.reserve(numPops);
  solvers_.reserve(numPops);

  std::vector<Eigen::Triplet<double>> coefficients(0);
  coefficients.reserve(sslib.getNumStats());

  // for each population (making this the outer loop seems to be the way to go)
  for(size_t i = 0; i < numPops; ++i)
  {
    std::string popId = sslib.asString(i);

    // for each stat in vector Y (going by rows of matrices_)
    for(auto it = std::begin(sslib->getStats()); it != std::end(sslib->getStats()); ++it)
    {
      std::string mom = (*it)->first; // full name of moment
      std::vector<std::string> splitMom = sslib.splitString(mom, "_"); // splits name by underscore

      size_t popIdCount = sslib.countInstances(splitMom[1], popId); // count of i in moment's name
      size_t row = it - std::begin(sslib->getStats()); // row index
      size_t col = 0; // column index

      if(splitMom[0] == "DD")
      {
        if(popIdCount == 2)
        {
          coefficients.push_back(Eigen::Triplet<double>(row, row -3.));

          col = sslib.indexLookup("Dz_" + popId + popId + popId);
          coefficients.push_back(Eigen::Triplet<double>(row, col, 1.));

          col = sslib.indexLookup("pi2_" + popId + popId + ";" + popId + popId);
          coefficients.push_back(Eigen::Triplet<double>(row, col, 1.));
        }

        else if(popIdCount == 1)
          coefficients.push_back(Eigen::Triplet<double>(row, row, -1.));
      }

      else if(splitMom[0] == "Dz")
      {
        if(splitMom[0] == popId) // if D_i_z**
        {
          if(popIdCount == 3)
          {
            coefficients.push_back(Eigen::Triplet<double>(row, row, -5.));

            col = sslib.indexLookup("DD_" + popId + popId);
            coefficients.push_back(Eigen::Triplet<double>(row, col, 4.));
          }

          else if(popIdCount == 2)
            coefficients.push_back(Eigen::Triplet<double>(row, row, -3.));

          else if(popIdCount == 1) // if D_i_z_xx where x != i
            coefficients.push_back(Eigen::Triplet<double>(row, row, -1.));
        }

        else // if D_x_z_** where x != i
        {
          if(popIdCount == 2)  // if D_x_z_ii where x != i
          {
            col = sslib.indexLookup("DD_" + popId + splitMom[1][0]);
            coefficients.push_back(Eigen::Triplet<double>(row, col, 4.));
          }
        }
      }

      else if(splitMom[0] == "pi2")
      {
        std::vector<std::string> splitPops = sslib.splitString(splitMom[1], ";"); // gets the 4 pop indices

        size_t countLeft = sslib.countInstances(splitPops[0], popId); // count of i before ';'
        size_t countRight = sslib.countInstances(splitPops[1], popId); // count of i after ';'

        if((countLeft + countRight) == 4)
        {
          coefficients.push_back(Eigen::Triplet<double>(row, row, -2.));

          col = sslib.indexLookup("Dz_" + popId + popId + popIds);
          coefficients.push_back(Eigen::Triplet<double>(row, col, 1.));
        }

        else if((countLeft == 2) || (countLeft == 2))
          coefficients.push_back(Eigen::Triplet<double>(row, row, -1.));
      }

      else if(splitMom[0] == "H")
        coefficients.push_back(Eigen::Triplet<double>(row, row, -1.));
    }

    Eigen::SparseMatrix<double> mat;
    mat.setFromTriplets(std::begin(coefficients), std::end(coefficients));
    mat.makeCompressed();
    matrices_.emplace_back(mat); // at the i-th position of vector, where i index the population

    Eigen::EigenSolver es(mat);
    solvers_.emplace_back(es); // TODO check if it's more efficient to store eigenvectors, eigenvalues, and eigenvectors ^ (-1)
  }
}

void Drift::updateMatrices_()
{
  std::string paramName = "";

  for(size_t i = 0; i < solvers_.size(); ++i)
  {
    paramName = "N_" + bpp::Textools::toString(i); // one per population

    double prevVal = prevParams_.getParameterValue(paramName); // old
    double newVal = getParameterValue(paramName); // new

    // we want something like this:
    solvers_[i].eigenvalues() *= (prevVal / newVal); // for Drift, it's inverted (relative to other operators) because we scale matrices by 1 / N
  }

  prevParams_.matchParametersValues(getParameters());
}
