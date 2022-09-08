/*
 * Authors: Gustavo V. Barroso
 * Created: 09/08/2022
 * Last modified: 08/09/2022
 *
 */

#include <Bpp/Text/TextTools.h>

#include "Drift.hpp"


void Drift::setUpMatrices_(const SumStatsLibrary& sslib)
{
  size_t numPops = sslib.getNumPops();

  matrices_.resize(numPops);
  eigenDec_.reserve(numPops);

  // for each population (making this the outer loop seems to be the way to go)
  for(size_t i = 0; i < numPops; ++i)
  {
    std::string popId = sslib.asString(i);

    matrices_[i] = Eigen::MatrixXd::Zero(sslib.getStats(), sslib.getStats()); // inits to 0 matrix

    // for each stat in vector Y (going by rows of matrices_)
    for(auto it = std::begin(sslib.getStats()); it != std::end(sslib.getStats()); ++it)
    {
      std::string mom = (*it)->first; // full name of moment
      std::vector<std::string> splitMom = sslib.splitString(mom, "_"); // splits name by underscore

      size_t popIdCount = sslib.countInstances(splitMom[1], popId); // count of i in moment's name
      size_t row = it - std::begin(sslib.getStats()); // row index
      size_t col = 0; // column index

      if(splitMom[0] == "DD")
      {
        if(popIdCount == 2)
        {
          matrices_[i](row, row) = -3.;

          col = sslib.indexLookup("Dz_" + popId + popId + popId);
          matrices_[i](row, col,) = 1.;

          col = sslib.indexLookup("pi2_" + popId + popId + ";" + popId + popId);
          matrices_[i](row, col) = 1.;
        }

        else if(popIdCount == 1)
          matrices_[i](row, row) = -1.;
      }

      else if(splitMom[0] == "Dz")
      {
        if(splitMom[0] == popId) // if D_i_z**
        {
          if(popIdCount == 3)
          {
            matrices_[i](row, row) = -5.;

            col = sslib.indexLookup("DD_" + popId + popId);
            matrices_[i](row, col) = 4.;
          }

          else if(popIdCount == 2)
            matrices_[i](row, row) = -3.;

          else if(popIdCount == 1) // if D_i_z_xx where x != i
            matrices_[i](row, row) = -1.;
        }

        else // if D_x_z_** where x != i
        {
          if(popIdCount == 2)  // if D_x_z_ii where x != i
          {
            col = sslib.indexLookup("DD_" + popId + splitMom[1][0]);
            matrices_[i](row, col) = 4.;
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
          matrices_[i](row, row) = -2.;

          col = sslib.indexLookup("Dz_" + popId + popId + popIds);
          matrices_[i](row, col) = 1.;
        }

        else if((countLeft == 2) || (countLeft == 2))
          matrices_[i](row, row) = -1.;
      }

      else if(splitMom[0] == "H")
        matrices_[i](row, row) = -1.;
    }

    // scale matrix by associated parameter value so that it
    matrices_[i] *= getParameterValue("N_" + bpp::TextTools::toString(i));

    EigenDecomposition ed(matrices_[i], exponent_);
    eigenDec_.emplace_back(ed);
  }
}

void Drift::updateMatrices_()
{
  std::string paramName = "";

  for(size_t i = 0; i < eigenDec_.size(); ++i) // one matrix / eigensolver per population (within each Epoch)
  {
    paramName = "N_" + bpp::Textools::toString(i);

    double prevVal = prevParams_.getParameterValue(paramName); // old
    double newVal = getParameterValue(paramName); // new
    double factor = std::pow(prevVal / newVal, exponent_); // for Drift, it's inverted (relative to other operators) because we scale matrices by 1 / N

    eigenDec_[i].setLambda(eigenDec_[i].lambda() * factor);
  }

  prevParams_.matchParametersValues(getParameters());
}
