/*
 * Authors: Gustavo V. Barroso
 * Created: 10/08/2022
 * Last modified: 08/09/2022
 *
 */


#include "Migration.hpp"


void Migration::setUpMatrices_(const SumStatsLibrary& sslib)
{
  size_t numPops = fetchNumPops();
  size_t index = 0; // matrix index, referring to the coefficients of the parameter m_ij, i != j

  matrices_.resize(numPops * (numPops - 1));
  eigenDec_.reserve(numPops * (numPops - 1));

  for(size_t i = 0; i < numPops; ++i)
  {
    std::string childPopId = bpp::TextTools::toString(i); // i in m_ij

    for(size_t j = 0; j < numPops; ++j)
    {
      if(i != j) // for each pair of populations (parameters m_ij, i != j)
      {
        matrices_[index] = Eigen::MatrixXd::Zero(sslib.getStats(), sslib.getStats()); // inits to 0 matrix
        std::string parentPopId = bpp::TextTools::toString(j); // j in m_ij

        // NOTE even though we deal with order_ >= 4 and have pi2(i,j;k,l) in stats,
        // we only need to loop over pops twice because the loop over stats takes care of the other pops (k,l)

        // for each stat in vector Y (rows of matrices_[index])
        for(auto it = std::begin(sslib->getStats()); it != std::end(sslib->getStats()); ++it)
        {
          std::string mom = *it->first; // full name of moment
          std::vector<std::string> splitMom = sslib.splitString(mom, "_"); // (see SumStatsLibrary::init)

          std::string p1, p2, p3, p4 = ""; // inits reference strings for population IDs

          size_t parentPopIdCount = sslib.countInstances(splitMom[1], parentPopId); // count of i in moment's name
          size_t row = it - std::begin(sslib->getStats()); // row index
          size_t col = 0; // column index

          // main diagonal initialization, this entry must exist in a Sparse matrix even if it remains zero because
          // it is important when converting from "delta" to "transition" matrix (see Operator::fetchCombinedMatrix()
          matrices_[index](row, row) = 0.;

          if(splitMom[0] == "DD")
          {
            matrices_[index](row, row) = -static_cast<double>(parentPopIdCount); // main diagonal

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
            matrices_[index](row, row) = -static_cast<double>(parentPopIdCount))); // main diagonal

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
                matrices_[index](row, col) = parentPopIdCount;
              }

              // if they both represent parent population
              else if(p2 == parentPopId && p3 == parentPopId)
              {
                // z_ij
                col = sslib.indexLookup("Dz_" + p1 + childPopId + parentPopId);
                matrices_[index](row, col) = 1.;

                // z_ji
                col = sslib.indexLookup("Dz_" + p1 + parentPopId + childPopId);
                matrices_[index](row, col) = 1.;
              }
            }

            // NOTE be mindful of interchangeable pi2 statistics (ordered vs non-ordered)
            else if(p1 == parentPopId) // if population index of D in Dz is j in m_ij
            {
              if(p2 == childPopId)
              {
                if(p3 == childPopId) // D_j_z_ii
                {
                  // the Dz cols
                  col = sslib.indexLookup("Dz_" + childPopId + childPopId + childPopId);
                  matrices_[index](row, col) = 1.;

                  // the pi2 cols
                  col = sslib.indexLookup("pi2_" + childPopId + childPopId + ";" + childPopId + childPopId); // append childPopId to left of p2 and p3
                  matrices_[index](row, col) = 4.;

                  col = sslib.indexLookup("pi2_" + childPopId + childPopId + ";" + childPopId + parentPopId); // switch right to parentPopId
                  matrices_[index](row, col) = -4.;

                  col = sslib.indexLookup("pi2_" + childPopId + parentPopId + ";" + childPopId + childPopId); // switch left to parentPopId
                  matrices_[index](row, col) = -4.;

                  col = sslib.indexLookup("pi2_" + childPopId + parentPopId + ";" + childPopId + parentPopId); // switch both
                  matrices_[index](row, col) = 4.;
                }

                else if(p3 == parentPopId) // D_j_z_ij
                {
                  // the Dz cols
                  col = sslib.indexLookup("Dz_" + parentPopId + childPopId + childPopId);
                  matrices_[index](row, col) = 1.;

                  col = sslib.indexLookup("Dz_" + childPopId + parentPopId + childPopId);
                  matrices_[index](row, col) = 1.;

                  // the pi2 cols
                  col = sslib.indexLookup("pi2_" + childPopId + childPopId + ";" + childPopId + parentPopId); // append childPopId to left of p2 and p3
                  matrices_[index](row, col) = 4.;

                  col = sslib.indexLookup("pi2_" + childPopId + childPopId + ";" + parentPopId + parentPopId); // switch right to parentPopId
                  matrices_[index](row, col) = -4.;

                  col = sslib.indexLookup("pi2_" + childPopId + parentPopId + ";" + childPopId + parentPopId); // switch left to parentPopId
                  matrices_[index](row, col) = -4.;

                  col = sslib.indexLookup("pi2_" + childPopId + parentPopId + ";" + parentPopId + parentPopId); // switch both
                  matrices_[index](row, col) = 4.;
                }
              }

              else if(p2 == parentPopId)
              {
                if(p3 == childPopId) // D_j_z_ji
                {
                  // the Dz cols
                  col = sslib.indexLookup("Dz_" + childPopId + parentPopId + childPopId);
                  matrices_[index](row, col) = 1.;

                  col = sslib.indexLookup("Dz_" + parentPopId + childPopId + childPopId);
                  matrices_[index](row, col) = 1.;

                  // the pi2 cols
                  col = sslib.indexLookup("pi2_" + childPopId + parentPopId + ";" + childPopId + childPopId); // append childPopId to left of p2 and p3
                  matrices_[index](row, col) = 4.;

                  col = sslib.indexLookup("pi2_" + childPopId + parentPopId + ";" + childPopId + parentPopId); // switch right to parentPopId
                  matrices_[index](row, col) = -4.;

                  col = sslib.indexLookup("pi2_" + parentPopId + parentPopId + ";" + childPopId + childPopId); // switch left to parentPopId
                  matrices_[index](row, col) = -4.;

                  col = sslib.indexLookup("pi2_" + parentPopId + parentPopId + ";" + childPopId + parentPopId); // switch both
                  matrices_[index](row, col) = 4.;
                }

                else if(p3 == parentPopId) // D_j_z_jj
                {
                  // the Dz cols
                  col = sslib.indexLookup("Dz_" + childPopId + parentPopId + parentPopId);
                  matrices_[index](row, col) = 1.;

                  col = sslib.indexLookup("Dz_" + parentPopId + childPopId + parentPopId);
                  matrices_[index](row, col) = 1.;

                  col = sslib.indexLookup("Dz_" + parentPopId + parentPopId + childPopId);
                  matrices_[index](row, col) = 1.;

                  // the pi2 cols
                  col = sslib.indexLookup("pi2_" + childPopId + parentPopId + ";" + childPopId + parentPopId); // append childPopId to left of p2 and p3
                  matrices_[index](row, col) = 4.;

                  col = sslib.indexLookup("pi2_" + childPopId + parentPopId + ";" + parentPopId + parentPopId); // switch right to parentPopId
                  matrices_[index](row, col) = -4.;

                  col = sslib.indexLookup("pi2_" + parentPopId + parentPopId + ";" + childPopId + parentPopId); // switch left to parentPopId
                  matrices_[index](row, col) = -4.;

                  col = sslib.indexLookup("pi2_" + parentPopId + parentPopId + ";" + parentPopId + parentPopId); // switch both
                  matrices_[index](row, col) = 4.;
                }
              }
            }
          }

          else if(splitMom[0] == "pi2")
          {
            std::vector<std::string> splitPops = sslib.splitString(splitMom[1], ";"); // splits name by semi-colon

            size_t countLeft = sslib.countInstances(splitPops[0], parentPopId); // count of j before ';'
            size_t countRight = sslib.countInstances(splitPops[1], parentPopId); // count of j after ';'

            // main diagonal init based on left pair
            matrices_[index](row, row) = -static_cast<double>(countLeft)));
            // updates based on right pair (duplicated elements are summed up by the Eigen method setFromTriplets)
            matrices_[index](row, row) = -static_cast<double>(countRight)));

            p1, p2, p3, p4 = childPopId; // resets reference strings

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
            matrices_[index](row, col) = countLeft;

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
            matrices_[index](row, col) = countRight;
          }

          else if(splitMom[0] == "H")
          {
            matrices_[index](row, row) = -static_cast<double>(parentPopIdCount); // main diagonal

            p1 = splitMom[1][0];
            p2 = splitMom[1][1];

            if(parentPopIdCount == 1)
            {
              if(p1 == childPopId || p2 == childPopId) // H_ij or H_ji
              {
                col = sslib.indexLookup("H_" + childPopId + childPopId);
                matrices_[index](row, col) = -static_cast<double>(parentPopIdCount);
              }
            }

            else if(parentPopIdCount == 2)
            {
              col = sslib.indexLookup("H_" + childPopId + parentPopId); // H_ij
              matrices_[index](row, col) = -static_cast<double>(parentPopIdCount); // -1 before compression

             /* H_ji -- before compression
              * col = sslib.indexLookup("H_" + parentPopId + childPopId);
              * matrices_[index](row, col) -= 1. * static_cast<double>(parentPopIdCount);
              */
            }
          }
        }

        matrices_[index] *= getParameterValue("m_" + bpp::TextTools::toString(i) + bpp::TextTools::toString(j));

        EigenDecomposition ed(matrices_[index], exponent);
        eigenDec_.emplace_back(ed);

        ++index;
      } // ends if(i != j)
    }
  }
}

void Migration::updateMatrices_()
{
  size_t numPops = fetchNumPops();

  size_t index = 0;
  std::string name = "";

  for(size_t i = 0; i < numPops; ++i)
  {
    for(size_t j = 0; j < numPops; ++j)
    {
      if(i != j)
      {
        name = "m_" + bpp::Textools::toString(i) + bpp::Textools::toString(j);

        double prevVal = prevParams_.getParameterValue(name);
        double newVal = getParameterValue(name);
        double factor = std::pow(newVal / prevVal, exponent_);

        eigenDec_[i].setLambda(eigenDec_[i].lambdaReal() * factor);

        ++index;
      }
    }
  }

  prevParams_.matchParametersValues(getParameters());
}

