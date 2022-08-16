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
  for(size_t i = 0; i < numPops; ++i)
  {
    std::string childPopId = asString(i); // i in m_ij

    for(size_t j = 0; j < numPops; ++j)
    {
      if(i != j) // for each pair of populations (parameters m_ij, i != j)
      {
        std::string parentPopId = asString(j); // j in m_ij

        // NOTE even though we deal with order_ = 4 and have pi2(i,j;k,l) in stats,
        // we only need to loop over pops twice because the loop over stats names takes care of the other pops (k,l)

        // for each stat in vector Y (going by rows of matrices_)
        for(auto it = std::begin(sslib->getStats()); it != std::end(sslib->getStats()); ++it)
        {
          std::string mom = *it->first; // full name of moment
          std::vector<std::string> splitMom = sslib.splitString(mom, "_"); // (cf SumStatsLibrary::init)

          std::string p1, p2, p3, p4 = childPopId; // starts reference strings for population IDs

          size_t parentPopIdCount = sslib.countInstances(splitMom[1], parentPopId); // count of i in moment's name
          size_t row = sslib.indexLookup(mom); // row index
          size_t col = 0; // column index

          if(splitMom[0] == "DD")
          {
            matrices_[i](row, row) = - static_cast<double>(parentPopIdCount); // main diagonal (-0.0 for 2nd-order children)

            if(parentPopIdCount == 1)
            {
              col = sslib.indexLookup("DD_" + p1 + p2);
              matrices_[i](row, col) = parentPopIdCount;
            }

            else if(parentPopIdCount == 2)
            {
              p2 = parentPopId;
              col = sslib.indexLookup("DD_" + p1 + p2); // WARNING there is no moment called DD_21 (only DD_12)
              matrices_[i](row, col) = parentPopIdCount;
            }
          }

          else if(splitMom[0] == "Dz")
          {

          }

          else if(splitMom[0] == "pi2")
          {
            std::vector<std::string> splitPops = sslib.splitString(mom, ";"); // splits name by semi-colon

            size_t countLeft = sslib.countInstances(splitPops[0], parentPopId); // count of j before ';'
            size_t countRight = sslib.countInstances(splitPops[1], parentPopId); // count of j after ';'

            // diagonal entry update based on left pair
            matrices_[i](row, row) = matrices_[i](row, row) - static_cast<double>(countLeft);
            // diagonal entry update based on right pair
            matrices_[i](row, row) = matrices_[i](row, row) - static_cast<double>(countRight);

            p1, p2, p3, p4 = childPopId; // re-starts reference strings

            // updates non-diagonal (positive) entries from left perspective
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

            // updates non-diagonal (positive) entries from right perspective
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
          }

          else if(splitMom[0] == "H")
          {
            matrices_[i](row, row) = -1. * static_cast<int>(parentPopIdCount); // diagonal entry

            sslib.indexLookup("H" + );
          }
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

