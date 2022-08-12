/*
 * Authors: Gustavo V. Barroso
 * Created: 09/08/2022
 * Last modified: 12/08/2022
 *
 */

#include <Bpp/Text/TextTools.h>

#include "Drift.hpp"


void Drift::setUpMatrix(const SumStatsLibrary& sslib)
{
  size_t numPops = sslib.getNumPops();

  // for each population (making this the outer loop seems to be the way to go)
  for(size_t i = 0; i < numPops; ++i)
  {
    std::string popId = asString(i);

    // for each stat in vector Y (going by rows of matrices_)
    for(auto it = std::begin(sslib->getStats()); it != std::end(sslib->getStats()); ++it)
    {
      std::string mom = *it->first; // full name of moment
      std::vector<std::string> splitMom = splitUnder_(mom); // splits name by underscore char

      size_t popIdCount = countInstances_(splitMom[1], popId); // count of i in moment's name
      size_t row = indexLookup_(mom); // row index
      size_t col = 0; // column index

      if(splitMom[0] == "DD")
      {
        if(popIdCount == 2)
        {
          matrices_[i](row, row) = -3.;

          col = indexLookup_("Dz_" + popId + popId + popId);
          matrices_[i](row, col) = 1.;

          col = indexLookup_("pi2_" + popId + popId + ";" + popId + popId);
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

            col = indexLookup_("DD_" + popId + popId);
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
            col = indexLookup_("DD_" + popId + splitMom[1][0]);
            matrices_[i](row, col) = 4.;
          }
        }
      }

      else if(splitMom[0] == "pi2")
      {
        std::vector<std::string> splitPops(0);
        boost::split(splitPops, splitMom[1], boost::is_any_of(";"));

        size_t popIdCountLeft = countInstances_(splitPops[0], popId); // count of i before ';'
        size_t popIdCountRight = countInstances_(splitPops[1], popId); // count of i after ';'

        if((popIdCountLeft + popIdCountRight) == 4)
        {
          matrices_[i](row, row) = -2.;

          col = indexLookup_("Dz_" + popId + popId + popIds);
          matrices_[i](row, col) = 1.;
        }

        else if((popIdCountLeft == 2) || (popIdCountLeft == 2))
          matrices_[i](row, row) = -1.;
      }

      else if(splitMom[0] == "H")
        matrices_[i](row, row) = -1.;
    }
  }
}

void Drift::update_()
{
  for(size_t i = 0; i < matrices_.size(); ++i)
  {
    std::string param = "N" + bpp::Textools::toString(i);

    matrices_[i] = (prevParams_[i] / getParameter(param)) * matrices_[i];
    prevParams_[i] = getParameter(param);
  }
}
