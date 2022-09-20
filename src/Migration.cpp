/*
 * Authors: Gustavo V. Barroso
 * Created: 10/08/2022
 * Last modified: 20/09/2022
 *
 */


#include "Migration.hpp"

void Migration::setUpMatrices_(const SumStatsLibrary& sslib)
{
  size_t numPops = fetchNumPops();
  size_t numStats = sslib.getNumStats();

  matrices_.reserve(numPops * (numPops - 1));

  for(size_t i = 0; i < numPops; ++i) // for i in m_ij (i)
  {
    for(size_t j = 0; j < numPops; ++j) // for j in m_ijn (j)
    {
      if(i != j) // if populations in pair are not the same
      {
        std::vector<Eigen::Triplet<double>> coefficients(0); // init sparse matrix coefficients
        coefficients.reserve(3 * numStats);

        // NOTE even though we deal with order_ >= 4 and have pi2(i,j;k,l) in stats,
        // we only need to loop over pops twice because the loop over stats takes care of the other pops (k,l)

        // for each stat in vector Y (rows of focal matrix)
        for(auto it = std::begin(sslib.getMoments()); it != std::end(sslib.getMoments()); ++it)
        {
          size_t p1, p2, p3, p4 = i; // inits reference indices for population IDs WARNING check if we want i

          size_t parentPopIdCount = it->countInstances(j); // count of j in moment's name WARNING check if we want i or j here
          size_t row = it - std::begin(sslib.getMoments()); // row index
          size_t col = 0; // column index

          if(it->getPrefix() == "DD")
          {
            coefficients.emplace_back(Eigen::Triplet<double>(row, row, -static_cast<double>(parentPopIdCount))); // main diagonal

            if(parentPopIdCount == 1)
            {
              col = sslib.findDdIndex(i, i); // WARNING check what we want for indices
              coefficients.emplace_back(Eigen::Triplet<double>(row, col, parentPopIdCount));
            }

            else if(parentPopIdCount == 2)
            {
              col = sslib.findDdIndex(i, j);
              coefficients.emplace_back(Eigen::Triplet<double>(row, col, parentPopIdCount));
            }
          }

          else if(it->getPrefix() == "Dz")
          {
            coefficients.emplace_back(Eigen::Triplet<double>(row, row, -static_cast<double>(parentPopIdCount))); // main diagonal

            p1 = it->getPopIndices()[0]; // D pop
            p2 = it->getPopIndices()[1]; // first z pop
            p3 = it->getPopIndices()[2]; // second z pop

            if(p1 == i) // if population index of D in Dz is i in m_ij
            {
              // if p2 and p3 represent child and parent populations, in any order
              if((p2 == i && p3 == j) || (p2 == j && p3 == i))
              {
                p2, p3 = i;

                col = sslib.findDzIndex(p1, p2, p3);
                coefficients.emplace_back(Eigen::Triplet<double>(row, col, parentPopIdCount));
              }

              // if they both represent parent population
              else if(p2 == j && p3 == j)
              {
                // z_ij
                col = sslib.findDzIndex(p1, i, j);
                coefficients.emplace_back(Eigen::Triplet<double>(row, col, 1.));

                // z_ji
                col = sslib.findDzIndex(p1, j, i);
                coefficients.emplace_back(Eigen::Triplet<double>(row, col, 1.));
              }
            }

            // NOTE be mindful of interchangeable pi2 statistics (ordered vs non-ordered)
            else if(p1 == j) // if population index of D in Dz is j in m_ij
            {
              if(p2 == i)
              {
                if(p3 == i) // D_j_z_ii
                {
                  // the Dz cols
                  col = sslib.findDzIndex(i, i, i);
                  coefficients.emplace_back(Eigen::Triplet<double>(row, col, 1.));

                  // the pi2 cols
                  col = sslib.findPi2Index(i, i, i, i); // append i to left of p2 and p3
                  coefficients.emplace_back(Eigen::Triplet<double>(row, col, 4.));

                  col = sslib.findPi2Index(i, i, i, j); // switch right to j
                  coefficients.emplace_back(Eigen::Triplet<double>(row, col, -4.));

                  col = sslib.findPi2Index(i, j, i, i); // switch left to j
                  coefficients.emplace_back(Eigen::Triplet<double>(row, col, -4.));

                  col = sslib.findPi2Index(i, j, i, j); // switch both
                  coefficients.emplace_back(Eigen::Triplet<double>(row, col, 4.));
                }

                else if(p3 == j) // D_j_z_ij
                {
                  // the Dz cols
                  col = sslib.findDzIndex(j, i, i);
                  coefficients.emplace_back(Eigen::Triplet<double>(row, col, 1.));

                  col = sslib.findDzIndex(i, j, i);
                  coefficients.emplace_back(Eigen::Triplet<double>(row, col, 1.));

                  // the pi2 cols
                  col = sslib.findPi2Index(i, i, i, j); // append i to left of p2 and p3
                  coefficients.emplace_back(Eigen::Triplet<double>(row, col, 4.));

                  col = sslib.findPi2Index(i, i, j, j); // switch right to j
                  coefficients.emplace_back(Eigen::Triplet<double>(row, col, -4.));

                  col = sslib.findPi2Index(i, j, i, j); // switch left to j
                  coefficients.emplace_back(Eigen::Triplet<double>(row, col, -4.));

                  col = sslib.findPi2Index(i, j, j, j); // switch both
                  coefficients.emplace_back(Eigen::Triplet<double>(row, col, 4.));
                }
              }

              else if(p2 == j)
              {
                if(p3 == i) // D_j_z_ji
                {
                  // the Dz cols
                  col = sslib.findDzIndex(i, j, i);
                  coefficients.emplace_back(Eigen::Triplet<double>(row, col, 1.));

                  col = sslib.findDzIndex(j, i, i);
                  coefficients.emplace_back(Eigen::Triplet<double>(row, col, 1.));

                  // the pi2 cols
                  col = sslib.findPi2Index(i, j, i, i); // append i to left of p2 and p3
                  coefficients.emplace_back(Eigen::Triplet<double>(row, col, 4.));

                  col = sslib.findPi2Index(i, j, i, j); // switch right to j
                  coefficients.emplace_back(Eigen::Triplet<double>(row, col, -4.));

                  col = sslib.findPi2Index(j, j, i, i); // switch left to j
                  coefficients.emplace_back(Eigen::Triplet<double>(row, col, -4.));

                  col = sslib.findPi2Index(j, j, i, j); // switch both
                  coefficients.emplace_back(Eigen::Triplet<double>(row, col, 4.));
                }

                else if(p3 == j) // D_j_z_jj
                {
                  // the Dz cols
                  col = sslib.findDzIndex(i, j, j);
                  coefficients.emplace_back(Eigen::Triplet<double>(row, col, 1.));

                  col = sslib.findDzIndex(j, i, j);
                  coefficients.emplace_back(Eigen::Triplet<double>(row, col, 1.));

                  col = sslib.findDzIndex(j, j, i);
                  coefficients.emplace_back(Eigen::Triplet<double>(row, col, 1.));

                  // the pi2 cols
                  col = sslib.findPi2Index(i, j, i, j); // append i to left of p2 and p3
                  coefficients.emplace_back(Eigen::Triplet<double>(row, col, 4.));

                  col = sslib.findPi2Index(i, j, j, j); // switch right to j
                  coefficients.emplace_back(Eigen::Triplet<double>(row, col, -4.));

                  col = sslib.findPi2Index(j, j, i, j); // switch left to j
                  coefficients.emplace_back(Eigen::Triplet<double>(row, col, -4.));

                  col = sslib.findPi2Index(j, j, j, j); // switch both
                  coefficients.emplace_back(Eigen::Triplet<double>(row, col, 4.));
                }
              }
            }
          }

          else if(it->getPrefix() == "pi2")
          {
            size_t countLeft = 0; // count of i before ';'
            size_t countRight = 0; // count of i after ';'

            if(it->getPopIndices()[0] == i)
             ++countLeft;

            if(it->getPopIndices()[1] == i)
             ++countLeft;

            if(it->getPopIndices()[2] == i)
              ++countRight;

            if(it->getPopIndices()[3] == i)
              ++countRight;

            // main diagonal init based on left pair
            coefficients.emplace_back(Eigen::Triplet<double>(row, row, -static_cast<double>(countLeft)));
            // updates based on right pair (duplicated elements are summed up by the Eigen method setFromTriplets)
            coefficients.emplace_back(Eigen::Triplet<double>(row, row, -static_cast<double>(countRight)));

            p1, p2, p3, p4 = i; // resets references

            // updates non-diagonal (positive) entries from left perspective
            if(countLeft == 1)
            {
              if(countRight == 1)
                p4 = j;

              else if(countRight == 2)
                p3, p4 = j;
            }

            else if(countLeft == 2)
            {
              p2 = j;

              if(countRight == 1)
                p4 = j;

              else if(countRight == 2)
                p3, p4 = j;
            }

            col = sslib.findPi2Index(p1, p2, p3, p4);
            coefficients.emplace_back(Eigen::Triplet<double>(row, col, countLeft));

            p1, p2, p3, p4 = i; // resets references

            // updates non-diagonal (positive) entries from right perspective
            if(countRight == 1)
            {
              if(countLeft == 1)
                p2 = j;

              else if(countLeft == 2)
                p1, p2 = j;
            }

            else if(countRight == 2)
            {
              p4 = j;

              if(countLeft == 1)
                p2 = j;

              else if(countLeft == 2)
                p1, p2 = j;
            }

            col = sslib.findPi2Index(p1, p2, p3, p4);
            coefficients.emplace_back(Eigen::Triplet<double>(row, col, countRight));
          }

          else if(it->getPrefix() == "H")
          {
            coefficients.emplace_back(Eigen::Triplet<double>(row, row, -static_cast<double>(parentPopIdCount))); // main diagonal

            p1 = it->getPopIndices()[0];
            p2 = it->getPopIndices()[1];

            if(parentPopIdCount == 1)
            {
              if(p1 == i || p2 == i) // H_ij or H_ji
              {
                col = sslib.findHetIndex(i, i);
                coefficients.emplace_back(Eigen::Triplet<double>(row, col, -static_cast<double>(parentPopIdCount)));
              }
            }

            else if(parentPopIdCount == 2)
            {
              col = sslib.findHetIndex(i, j); // H_ij
              coefficients.emplace_back(Eigen::Triplet<double>(row, col, -static_cast<double>(parentPopIdCount))); // -1 before compression

             /* H_ji -- before compression
              * col = sslib.indexLookup("H_" + j + i);
              * coefficients.emplace_back(Eigen::Triplet<double>(row, col) -= 1. * static_cast<double>(parentPopIdCount);
              */
            }
          }
        }

        Eigen::SparseMatrix<double> mat(numStats, numStats);
        mat.setFromTriplets(std::begin(coefficients), std::end(coefficients));
        mat.makeCompressed();
        mat *= getParameterValue("m_" + bpp::TextTools::toString(i) + bpp::TextTools::toString(j)); // scales because of updateMatrices_()
        matrices_.emplace_back(mat);
      } // ends if(i != j)
    }
  }
}

void Migration::updateMatrices_()
{
  size_t numPops = fetchNumPops();

  size_t index = 0;
  std::string paramName = "";

  for(size_t i = 0; i < numPops; ++i)
  {
    for(size_t j = 0; j < numPops; ++j)
    {
      if(i != j)
      {
        paramName = "m_" + bpp::TextTools::toString(i) + bpp::TextTools::toString(j);

        double prevVal = prevParams_.getParameterValue(paramName);
        double newVal = getParameterValue(paramName);

        matrices_[index] *= (newVal / prevVal);

        ++index;
      }
    }
  }

  prevParams_.matchParametersValues(getParameters());
}

