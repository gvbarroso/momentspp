/*
 * Authors: Gustavo V. Barroso
 * Created: 10/08/2022
 * Last modified: 23/09/2022
 *
 */


#include "Migration.hpp"

void Migration::setUpMatrices_(const SumStatsLibrary& sslib)
{
  size_t numPops = fetchNumPops();
  size_t numStats = sslib.getNumStats();

  matrices_.reserve(numPops * (numPops - 1));

  for(size_t i = 0; i < numPops; ++i) // for i in m_ij (children pop id)
  {
    for(size_t j = 0; j < numPops; ++j) // for j in m_ij (parent pop id)
    {
      if(i != j) // if populations in pair are not the same, fills in focal matrix (scaled by m_ij)
      {
        std::vector<Eigen::Triplet<double>> coefficients(0); // init sparse matrix coefficients
        coefficients.reserve(3 * numStats);

        // NOTE even though we have pi2(i,j;k,l) in stats,
        // we only need to loop over pops twice because the loop over stats takes care of the other pops (k,l)

        // for each stat in vector Y (rows of focal matrix), lexicographically sorted based on Moments names
        for(auto it = std::begin(sslib.getMoments()); it != std::end(sslib.getMoments()); ++it)
        {
          size_t parentPopIdCount = it->countInstances(j); // count of j in Moment's name
          size_t row = it - std::begin(sslib.getMoments()); // row index
          size_t col = 0; // column index

          if(it->getPrefix() == "DD")
          {
            coefficients.emplace_back(Eigen::Triplet<double>(row, row, -parentPopIdCount)); // main diagonal

            if(parentPopIdCount == 1)
            {
              col = sslib.findDdIndex(i, i);
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
            coefficients.emplace_back(Eigen::Triplet<double>(row, row, -parentPopIdCount)); // main diagonal

            size_t p1 = it->getPopIndices()[0]; // D pop
            size_t p2 = it->getPopIndices()[1]; // first z pop
            size_t p3 = it->getPopIndices()[2]; // second z pop

            if(p1 == i) // if population index of D in Dz*** is i in m_ij
            {
              // if p2 and p3 represent child and parent populations, in any order
              if((p2 == i && p3 == j) || (p2 == j && p3 == i))
              {
                col = sslib.findDzIndex(i, i, i);
                coefficients.emplace_back(Eigen::Triplet<double>(row, col, parentPopIdCount));
              }

              // if they both represent parent population
              else if(p2 == j && p3 == j)
              {
                // z_ij
                col = sslib.findDzIndex(i, i, j);
                coefficients.emplace_back(Eigen::Triplet<double>(row, col, 1.));

                // z_ji
                col = sslib.findDzIndex(i, j, i);
                coefficients.emplace_back(Eigen::Triplet<double>(row, col, 1.));
              }
            }

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
                  // imagine starting with pop indices p2 and p3 on each side of ';' character in pi2(**;**)
                  // now append i to the right of both p2 and p3

                  col = sslib.findPi2Index(i, i, i, i); // as is
                  coefficients.emplace_back(Eigen::Triplet<double>(row, col, 4.));

                  col = sslib.findPi2Index(i, i, i, j); // switch right appendix to j
                  coefficients.emplace_back(Eigen::Triplet<double>(row, col, -4.));

                  col = sslib.findPi2Index(i, j, i, i); // switch left appendix to j
                  coefficients.emplace_back(Eigen::Triplet<double>(row, col, -4.));

                  col = sslib.findPi2Index(i, j, i, j); // switch both
                  coefficients.emplace_back(Eigen::Triplet<double>(row, col, 4.));
                }

                else if(p3 == j) // D_j_z_ij
                {
                  // the Dz cols
                  col = sslib.findDzIndex(j, i, i);
                  coefficients.emplace_back(Eigen::Triplet<double>(row, col, 1.));

                  col = sslib.findDzIndex(i, i, j);
                  coefficients.emplace_back(Eigen::Triplet<double>(row, col, 1.));

                  // the pi2 cols
                  // imagine starting with pop indices p2 and p3 on each side of ';' character in pi2(**;**)
                  // now append i to right the of p2 and to the left of p3

                  col = sslib.findPi2Index(i, i, i, j); // as is
                  coefficients.emplace_back(Eigen::Triplet<double>(row, col, 4.));

                  col = sslib.findPi2Index(i, i, j, j); // switch right appendix to j
                  coefficients.emplace_back(Eigen::Triplet<double>(row, col, -4.));

                  col = sslib.findPi2Index(i, j, i, j); // switch left appendix to j
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
                  // imagine starting with pop indices p2 and p3 on each side of ';' character in pi2(**;**)
                  // now append i to left of p2 and right of p3

                  col = sslib.findPi2Index(i, j, i, i); // as is
                  coefficients.emplace_back(Eigen::Triplet<double>(row, col, 4.));

                  col = sslib.findPi2Index(i, j, i, j); // switch right appendix to j
                  coefficients.emplace_back(Eigen::Triplet<double>(row, col, -4.));

                  col = sslib.findPi2Index(j, j, i, i); // switch left appendix to j
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
                  // imagine starting with pop indices p2 and p3 on each side of ';' character in pi2(**;**)
                  // now append i to left of both p2 and p3

                  col = sslib.findPi2Index(i, j, i, j); // as is
                  coefficients.emplace_back(Eigen::Triplet<double>(row, col, 4.));

                  col = sslib.findPi2Index(i, j, j, j); // switch right appendix to j
                  coefficients.emplace_back(Eigen::Triplet<double>(row, col, -4.));

                  col = sslib.findPi2Index(j, j, i, j); // switch left appendix to j
                  coefficients.emplace_back(Eigen::Triplet<double>(row, col, -4.));

                  col = sslib.findPi2Index(j, j, j, j); // switch both
                  coefficients.emplace_back(Eigen::Triplet<double>(row, col, 4.));
                }
              }
            }
          }

          else if(it->getPrefix() == "pi2")
          {
            size_t countLeft = 0; // count of j before ';' character in pi2(**;**)
            size_t countRight = 0; // count of j after ';' character in pi2(**;**)

            if(it->getPopIndices()[0] == j)
             ++countLeft;

            if(it->getPopIndices()[1] == j)
             ++countLeft;

            if(it->getPopIndices()[2] == j)
              ++countRight;

            if(it->getPopIndices()[3] == j)
              ++countRight;

            // main diagonal
            coefficients.emplace_back(Eigen::Triplet<double>(row, row, -countLeft - countRight));

            // updates off main-diagonal entries from left perspective
            size_t p1 = i;
            size_t p2 = i;
            size_t p3 = i;
            size_t p4 = i;

            if(countLeft == 1)
            {
              if(countRight == 1)
                p4 = j;

              else if(countRight == 2)
              {
                p3 = j;
                p4 = j;
              }
            }

            else if(countLeft == 2)
            {
              p2 = j;

              if(countRight == 1)
                p4 = j;

              else if(countRight == 2)
              {
                p3 = j;
                p4 = j;
              }
            }

            col = sslib.findPi2Index(p1, p2, p3, p4);
            coefficients.emplace_back(Eigen::Triplet<double>(row, col, countLeft));

            // updates off main-diagonal entries from right perspective
            p1 = i;
            p2 = i;
            p3 = i;
            p4 = i;

            if(countRight == 1)
            {
              if(countLeft == 1)
                p2 = j;

              else if(countLeft == 2)
              {
                p1 = j;
                p2 = j;
              }
            }

            else if(countRight == 2)
            {
              p4 = j;

              if(countLeft == 1)
                p2 = j;

              else if(countLeft == 2)
              {
                p1 = j;
                p2 = j;
              }
            }

            col = sslib.findPi2Index(p1, p2, p3, p4);
            coefficients.emplace_back(Eigen::Triplet<double>(row, col, countRight));
          }

          else if(it->getPrefix() == "H")
          {
            coefficients.emplace_back(Eigen::Triplet<double>(row, row, -parentPopIdCount)); // main diagonal

            size_t p1 = it->getPopIndices()[0];
            size_t p2 = it->getPopIndices()[1];

            if(parentPopIdCount == 1)
            {
              if(p1 == i || p2 == i) // H_ij or H_ji
              {
                col = sslib.findHetIndex(i, i);
                coefficients.emplace_back(Eigen::Triplet<double>(row, col, parentPopIdCount));
              }
            }

            else if(parentPopIdCount == 2)
            {
              col = sslib.findHetIndex(i, j);
              coefficients.emplace_back(Eigen::Triplet<double>(row, col, -1.));

              col = sslib.findHetIndex(j, i);
              coefficients.emplace_back(Eigen::Triplet<double>(row, col, -1.));
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

