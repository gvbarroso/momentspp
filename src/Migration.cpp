/*
 * Authors: Gustavo V. Barroso
 * Created: 10/08/2022
 * Last modified: 18/08/2022
 *
 */


#include "Migration.hpp"


void Migration::setUpMatrix(const SumStatsLibrary& sslib)
{
  size_t numPops = sslib.getNumPops();
  size_t index = 0; // matrix index, referring to the coefficients of the parameter m_ij, i != j

  // filling in focal matrix matrices_[index]
  for(size_t i = 0; i < numPops; ++i)
  {
    std::string childPopId = sslib.asString(i); // i in m_ij

    for(size_t j = 0; j < numPops; ++j)
    {
      if(i != j) // for each pair of populations (parameters m_ij, i != j)
      {
        std::string parentPopId = sslib.asString(j); // j in m_ij

        // NOTE even though we deal with order_ = 4 and have pi2(i,j;k,l) in stats,
        // we only need to loop over pops twice because the loop over stats' names takes care of the other pops (k,l)

        // for each stat in vector Y (rows of matrices_[index])
        for(auto it = std::begin(sslib->getStats()); it != std::end(sslib->getStats()); ++it)
        {
          // TODO add another loop over vector Y to represent cols of matrices_[index] to generalize to multiple populations?

          std::string mom = *it->first; // full name of moment
          std::vector<std::string> splitMom = sslib.splitString(mom, "_"); // (see SumStatsLibrary::init)

          std::string p1, p2, p3, p4 = childPopId; // inits reference strings for population IDs

          size_t parentPopIdCount = sslib.countInstances(splitMom[1], parentPopId); // count of i in moment's name
          size_t row = sslib.indexLookup(mom); // row index
          size_t col = 0; // column index

          if(splitMom[0] == "DD")
          {
            matrices_[index](row, row) -= static_cast<double>(parentPopIdCount); // main diagonal

            if(parentPopIdCount == 1)
            {
              col = sslib.indexLookup("DD_" + p1 + p2);
              matrices_[index](row, col) = parentPopIdCount;
            }

            else if(parentPopIdCount == 2)
            {
              p2 = parentPopId;
              col = sslib.indexLookup("DD_" + p1 + p2);
              matrices_[index](row, col) = parentPopIdCount;
            }
          }

          else if(splitMom[0] == "Dz")
          {
            matrices_[index](row, row) -= static_cast<double>(parentPopIdCount); // main diagonal

            p1 = splitMom[1][0]; // D pop
            p2 = splitMom[1][1]; // first z pop
            p3 = splitMom[1][2]; // second z pop

            if(p1 == childPopId) // if population index of D in Dz is i in m_ij
            {
              // if p2 and p3 represent child and parent populations, in any order
              if((p2 == childPopId && p3 == parentPopId) || (p2 == parentPopId && p3 == childPopId))
              {
                p2, p3 = childPopId;

                col = sslib.indexLookup("Dz_" + p1 + p2 + p3);
                matrices_[index](row, col) += parentPopIdCount;
              }

              // if they both represent parent population
              else if(p2 == parentPopId && p3 == parentPopId)
              {
                // z_ij
                col = sslib.indexLookup("Dz_" + p1 + childPopId + parentPopId);
                matrices_[index](row, col) += 1.;

                // z_ji
                col = sslib.indexLookup("Dz_" + p1 + parentPopId + childPopId);
                matrices_[index](row, col) += 1.;
              }
            }

            else if(p1 == parentPopId) // if population index of D in Dz is j in m_ij
            {
              // helper variables for pi2 entries
              std::string pair1 = childPopId + p2;
              std::string pair2 = childPopId + p3;

              if(p2 == childPopId)
              {
                if(p3 == childPopId) // D_j_z_ii
                {
                  // the Dz cols
                  col = sslib.indexLookup("Dz_" + childPopId + childPopId + childPopId);
                  matrices_[index](row, col) += 1.;

                  // the pi2 cols
                  col = sslib.indexLookup("pi2_" + pair1 + ";" + pair2);
                  matrices_[index](row, col) += 4.;

                }

                else if(p3 == parentPopId) // D_j_z_ij
                {
                  // the Dz cols
                  col = sslib.indexLookup("Dz_" + parentPopId + childPopId + childPopId);
                  matrices_[index](row, col) += 1.;

                  col = sslib.indexLookup("Dz_" + childPopId + parentPopId + childPopId);
                  matrices_[index](row, col) += 1.;

                  // the pi2 cols
                  // ...
                }
              }

              else if(p2 == parentPopId)
              {
                if(p3 == childPopId) // D_j_z_ji
                {
                  // the Dz cols
                  col = sslib.indexLookup("Dz_" + childPopId + parentPopId + childPopId);
                  matrices_[index](row, col) += 1.;

                  col = sslib.indexLookup("Dz_" + parentPopId + childPopId + childPopId);
                  matrices_[index](row, col) += 1.;

                  // the pi2 cols
                  // ...
                }

                else if(p3 == parentPopId) // D_j_z_jj
                {
                  // the Dz cols
                  col = sslib.indexLookup("Dz_" + childPopId + parentPopId + parentPopId);
                  matrices_[index](row, col) += 1.;

                  col = sslib.indexLookup("Dz_" + parentPopId + childPopId + parentPopId);
                  matrices_[index](row, col) += 1.;

                  col = sslib.indexLookup("Dz_" + parentPopId + parentPopId + childPopId);
                  matrices_[index](row, col) += 1.;

                  // the pi2 cols
                  // ...
                }
              }
            }
          }

          else if(splitMom[0] == "pi2")
          {
            std::vector<std::string> splitPops = sslib.splitString(mom, ";"); // splits name by semi-colon

            size_t countLeft = sslib.countInstances(splitPops[0], parentPopId); // count of j before ';'
            size_t countRight = sslib.countInstances(splitPops[1], parentPopId); // count of j after ';'

            // main diagonal entry update based on left pair
            matrices_[index](row, row) -= static_cast<double>(countLeft);
            // main diagonal entry update based on right pair
            matrices_[index](row, row) -= static_cast<double>(countRight);

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
            matrices_[index](row, col) = matrices_[index](row, col) + countLeft;

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
            matrices_[index](row, col) = matrices_[index](row, col) + countRight;
          }

          else if(splitMom[0] == "H")
          {
            matrices_[index](row, row) = -1. * static_cast<int>(parentPopIdCount); // diagonal entry

            sslib.indexLookup("H" + );
          }
        }

        ++index;
      } // end if(i != j)
    }
  }
}

void Migration::update_()
{
  // this is a weird-looking but fun way to get the number of populations P from the raw value of P choose 2 ( == matrices_.size())
  int numPops = 0; // we want the positive solution of the quadratic equation P^2 - P - matrices_.size() = 0
  int binCoeff = static_cast<int>(matrices_.size()); // raw value of P choose 2
  for(int i = 2; i < binCoeff; ++i) // we never hit the upper bound of this loop but whatever
  {
    if(i * (1 - i) == -binCoeff)  // guaranteed to find if matrices_.size() was built correctly
    {
      numPops = i;
      break;
    }
  }

  size_t index = 0;
  std::string paramName = "";

  for(size_t i = 0; i < numPops; ++i)
  {
    for(size_t j = 0; j < numPops; ++j)
    {
      if(i != j)
      {
        paramName = "m_" + bpp::Textools::toString(i) + bpp::Textools::toString(j);

        double prevVal = prevParams_.getParameterValue(paramName);
        double newVal = getParameterValue(paramName); // from within itself

        matrices_[index] *= (newVal / prevVal);

        prevParams_.setParameterValue(paramName, newVal); // moves along

        ++index;
      }
    }
  }

  combineMatrices_();
}

