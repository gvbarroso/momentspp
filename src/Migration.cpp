/*
 * Authors: Gustavo V. Barroso
 * Created: 10/08/2022
 * Last modified: 16/08/2022
 *
 */


#include "Migration.hpp"


void Migration::setUpMatrix(const SumStatsLibrary& sslib)
{
  size_t numPops = sslib.getNumPops();

  // filling in focal matrix matrices_[i]
  for(size_t i = 0; i < numPops - 1; ++i)
  {
    std::string childPopId = asString(i); // i in m_ij

    // for each pair of populations (parameters m_ij, i != j)
    for(size_t j = i + 1; j < numPops; ++j)
    {
      std::string parentPopId = asString(j); // j in m_ij

      // NOTE even though we deal with order_ = 4 and have pi2(i,j;k,l) in stats,
      // we only need to loop over pops twice because the loop over stats names takes care of the other pops (k,l)

      // for each stat in vector Y (going by rows of matrices_)
      for(auto it = std::begin(sslib->getStats()); it != std::end(sslib->getStats()); ++it)
      {
        std::string mom = *it->first; // full name of moment
        std::vector<std::string> splitMom = sslib.splitString(mom, "_"); // (cf SumStatsLibrary::init)

        size_t parentPopIdCount = sslib.countInstances(splitMom[1], parentPopId); // count of i in moment's name
        size_t row = sslib.indexLookup(mom); // row index
        size_t col = 0; // column index

        if(splitMom[0] == "DD")
        {

        }

        else if(splitMom[0] == "Dz")
        {

        }

        else if(splitMom[0] == "pi2")
        {
          std::vector<std::string> splitPops = sslib.splitString(mom, ";"); // splits name by semi-colon

          size_t countLeft = sslib.countInstances(splitPops[0], parentPopId); // count of i before ';'
          size_t countRight = sslib.countInstances(splitPops[1], parentPopId); // count of i after ';'

          // diagonal entry update based on left pair
          matrices_[i](row, row) = matrices_[i](row, row) - static_cast<double>(countLeft);
          // diagonal entry update based on right pair
          matrices_[i](row, row) = matrices_[i](row, row) - static_cast<double>(countRight);

          std::string p1, p2, p3, p4 = childPopId; // starts reference strings

          // left perspective
          if(countLeft == 1)
          {
            if(countRight == 1)
              p4 = parentPopId;

            else if(countRight == 2)
              p3, p4 = parentPopId;
          }

          else if(countLeft == 2)
          {
            p2 = parentPopId;

            if(countRight == 1)
              p4 = parentPopId;

            else if(countRight == 2)
              p3, p4 = parentPopId;
          }

          col = sslib.indexLookup("pi2_" + p1 + p2 + ";" + p3 + p4);
          matrices_[i](row, col) = matrices_[i](row, col) + countLeft;

          p1, p2, p3, p4 = childPopId; // resets reference strings

          // right perspective
          if(countRight == 1)
          {
            if(countLeft == 1)
              p2 = parentPopId;

            else if(countLeft == 2)
              p1, p2 = parentPopId;
          }

          else if(countRight == 2)
          {
            p4 = parentPopId;

            if(countLeft == 1)
              p2 = parentPopId;

            else if(countLeft == 2)
              p1, p2 = parentPopId;
          }

          col = sslib.indexLookup("pi2_" + p1 + p2 + ";" + p3 + p4);
          matrices_[i](row, col) = matrices_[i](row, col) + countRight;


          /* non-diagonal update (positive entries) => we need control flow over the counts to look up the correct indices
          // left perspective
          if(countLeft == 1)
          {
            if(countRight == 0)
            {
              col = sslib.indexLookup("pi2_" + childPopId + childPopId + ";" + childPopId + childPopId);
              matrices_[i](row, col) = matrices_[i](row, col) + countLeft;
            }

            else if(countRight == 1)
            {
              col = sslib.indexLookup("pi2_" + childPopId + childPopId + ";" + childPopId + parentPopId);
              matrices_[i](row, col) = matrices_[i](row, col) + countLeft;
            }

            else if(countRight == 2)
            {
              col = sslib.indexLookup("pi2_" + childPopId + childPopId + ";" + parentPopId + parentPopId);
              matrices_[i](row, col) = matrices_[i](row, col) + countLeft;
            }
          }

          else if(countLeft == 2)
          {
            if(countRight == 0)
            {
              col = sslib.indexLookup("pi2_" + childPopId + parentPopId + ";" + childPopId + childPopId);
              matrices_[i](row, col) = matrices_[i](row, col) + countLeft;
            }

            else if(countRight == 1)
            {
              col = sslib.indexLookup("pi2_" + childPopId + parentPopId + ";" + childPopId + parentPopId);
              matrices_[i](row, col) = matrices_[i](row, col) + countLeft;
            }

            else if(countRight == 2)
            {
              col = sslib.indexLookup("pi2_" + childPopId + parentPopId + ";" + parentPopId + parentPopId);
              matrices_[i](row, col) = matrices_[i](row, col) + countLeft;
            }
          }

          // right perspective
          if(countRight == 1)
          {
            if(countLeft == 0)
            {
              col = sslib.indexLookup("pi2_" + childPopId + childPopId + ";" + childPopId + childPopId);
              matrices_[i](row, col) = matrices_[i](row, col) + countRight;
            }

            else if(countLeft == 1)
            {
              col = sslib.indexLookup("pi2_" + childPopId + parentPopId + ";" + childPopId + childPopId);
              matrices_[i](row, col) = matrices_[i](row, col) + countRight;
            }

            else if(countLeft == 2)
            {
              col = sslib.indexLookup("pi2_" + parentPopId + parentPopId + ";" + childPopId + childPopId);
              matrices_[i](row, col) = matrices_[i](row, col) + countRight;
            }
          }

          else if(countRight == 2)
          {
            if(countLeft == 0)
            {
              col = sslib.indexLookup("pi2_" + childPopId + childPopId + ";" + childPopId + parentPopId);
              matrices_[i](row, col) = matrices_[i](row, col) + countRight;
            }

            else if(countLeft == 1)
            {
              col = sslib.indexLookup("pi2_" + childPopId + parentPopId + ";" + childPopId + parentPopId);
              matrices_[i](row, col) = matrices_[i](row, col) + countRight;
            }

            else if(countLeft == 2)
            {
              col = sslib.indexLookup("pi2_" + parentPopId + parentPopId + ";" + childPopId + parentPopId);
              matrices_[i](row, col) = matrices_[i](row, col) + countRight;
            }
          }
        }*/

        else if(splitMom[0] == "H")
        {
          matrices_[i](row, row) = -1. * static_cast<int>(parentPopIdCount); // diagonal entry

          sslib.indexLookup("H" + );
        }
      }
    }
  }
}

void Migration::update_()
{
  matrix_ =  (getParameter() / prevParam_) * matrix_;
  prevParam_ = getParameter();
}

