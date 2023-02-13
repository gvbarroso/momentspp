/*
 * Authors: Gustavo V. Barroso
 * Created: 10/08/2022
 * Last modified: 08/02/2023
 *
 */


#include "Migration.hpp"


void Migration::setUpMatrices_(const SumStatsLibrary& sslib)
{
  // m_ij is the forward migration rate from pop i to pop j (backwards, the prob that lineage in j has parent in i)
  size_t numPops = fetchNumPops();
  size_t numStats = sslib.getNumStats();

  matrices_.reserve(numPops * (numPops - 1));

  for(size_t i = 0; i < numPops; ++i) // for i in m_ij
  {
    for(size_t j = 0; j < numPops; ++j) // for j in m_ij
    {
      if(i != j) // if populations in pair are not the same, fills in coefficients in focal matrix
      {
        std::vector<Eigen::Triplet<double>> coeffs(0); // init sparse matrix coeffs
        coeffs.reserve(3 * numStats);

        // although we pi2(i,j;k,l) moments, we need only to loop over pops twice because the loop over stats takes care of the k,l pops
        // for each stat in vector Y (rows of focal matrix), lexicographically sorted based on Moments names
        for(auto it = std::begin(sslib.getMoments()); it != std::end(sslib.getMoments()); ++it)
        {
          //(*it)->printAttributes(std::cout); // NOTE DEBUG

          int row = it - std::begin(sslib.getMoments()); // row index
          int col = -1; // inits column index to out-of-bounds
          int childPopIdCount = static_cast<int>((*it)->countInstances(j));
          if((*it)->getPrefix() == "DD")
          {
            if(childPopIdCount != 0) // not to populate the sparse matrix with unnecessary zeros
              coeffs.emplace_back(Eigen::Triplet<double>(row, row, -childPopIdCount));

            if(childPopIdCount == 1)
            {
              col = sslib.findDdIndex(i, i);
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));
            }

            else if(childPopIdCount == 2)
            {
              // look for covariance in D w.r.t to parent population (index i)
              col = sslib.findDdIndex(i, j);
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));

              col = sslib.findDdIndex(j, i);
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));
            }
          }

          else if((*it)->getPrefix() == "Dz")
          {
            if(childPopIdCount != 0) // not to populate the sparse matrix with unnecessary zeros
              coeffs.emplace_back(Eigen::Triplet<double>(row, row, -childPopIdCount));

            size_t p1 = (*it)->getPopIndices()[0]; // D pop
            size_t p2 = (*it)->getPopIndices()[1]; // first z pop
            size_t p3 = (*it)->getPopIndices()[2]; // second z pop

            if(p1 == i) // D_i_z_**
            {
              // D_i_z_ij || D_i_z_ji
              if((p2 == i && p3 == j) || (p2 == j && p3 == i))
              {
                col = sslib.findDzIndex(i, i, i);
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));
              }

              // D_i_z_jj
              else if(p2 == j && p3 == j)
              {
                col = sslib.findDzIndex(i, i, j);
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));

                col = sslib.findDzIndex(i, j, i);
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));
              }
            }

            else if(p1 == j) // D_j_z_**
            {
              if(p2 == i)
              {
                if(p3 == i) // D_j_z_ii
                {
                  // the Dz cols
                  col = sslib.findDzIndex(i, i, i);
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));

                  // the pi2 cols (minding permutations of derived/ancestral encoded in the suffixes, see Pi2Moment class)
                  // imagine starting with pop indices p2 and p3 on each side of ';' character in pi2(**;**) -> [i*;*i]
                  // now append i to the right of both p2 and p3

                  // as is
                  col = sslib.findPi2Index(i, i, i, i, "A");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));

                  col = sslib.findPi2Index(i, i, i, i, "B");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));

                  col = sslib.findPi2Index(i, i, i, i, "C");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));

                  col = sslib.findPi2Index(i, i, i, i, "D");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));

                  // switch right appendix to j
                  col = sslib.findPi2Index(i, i, i, j, "A");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./2.));

                  col = sslib.findPi2Index(i, i, i, j, "B");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./2.));

                  col = sslib.findPi2Index(i, i, i, j, "C");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./2.));

                  col = sslib.findPi2Index(i, i, i, j, "D");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./2.));

                  // permute pop in right locus
                  col = sslib.findPi2Index(i, i, j, i, "A");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./2.));

                  col = sslib.findPi2Index(i, i, j, i, "B");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./2.));

                  col = sslib.findPi2Index(i, i, j, i, "C");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./2.));

                  col = sslib.findPi2Index(i, i, j, i, "D");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./2.));

                  // switch left appendix to j
                  col = sslib.findPi2Index(i, j, i, i, "A");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./2.));

                  col = sslib.findPi2Index(i, j, i, i, "B");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./2.));

                  col = sslib.findPi2Index(i, j, i, i, "C");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./2.));

                  col = sslib.findPi2Index(i, j, i, i, "D");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./2.));

                  // permute pop in left locus
                  col = sslib.findPi2Index(j, i, i, i, "A");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./2.));

                  col = sslib.findPi2Index(j, i, i, i, "B");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./2.));

                  col = sslib.findPi2Index(j, i, i, i, "C");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./2.));

                  col = sslib.findPi2Index(j, i, i, i, "D");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./2.));

                  // switch both
                  col = sslib.findPi2Index(i, j, i, j, "A");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./4.));

                  col = sslib.findPi2Index(i, j, i, j, "B");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./4.));

                  col = sslib.findPi2Index(i, j, i, j, "C");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./4.));

                  col = sslib.findPi2Index(i, j, i, j, "D");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./4.));

                  // permute pop in right locus
                  col = sslib.findPi2Index(i, j, j, i, "A");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./4.));

                  col = sslib.findPi2Index(i, j, j, i, "B");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./4.));

                  col = sslib.findPi2Index(i, j, j, i, "C");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./4.));

                  col = sslib.findPi2Index(i, j, j, i, "D");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./4.));

                  // permute pop in left locus
                  col = sslib.findPi2Index(j, i, i, j, "A");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./4.));

                  col = sslib.findPi2Index(j, i, i, j, "B");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./4.));

                  col = sslib.findPi2Index(j, i, i, j, "C");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./4.));

                  col = sslib.findPi2Index(j, i, i, j, "D");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./4.));

                  // permute pops in both loci
                  col = sslib.findPi2Index(j, i, j, i, "A");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./4.));

                  col = sslib.findPi2Index(j, i, j, i, "B");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./4.));

                  col = sslib.findPi2Index(j, i, j, i, "C");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./4.));

                  col = sslib.findPi2Index(j, i, j, i, "D");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./4.));
                }

                else if(p3 == j) // D_j_z_ij
                {
                  // the Dz cols
                  col = sslib.findDzIndex(j, i, i);
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));

                  col = sslib.findDzIndex(i, i, j);
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./2.));

                  col = sslib.findDzIndex(i, j, i);
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./2.));

                  // the pi2 cols
                  // imagine starting with pop indices p2 and p3 on each side of ';' character in pi2(**;**) -> [i*;*j]
                  // now append i to right the of p2 and to the left of p3

                  // as is
                  col = sslib.findPi2Index(i, i, i, j, "A");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./2.));

                  col = sslib.findPi2Index(i, i, i, j, "B");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./2.));

                  col = sslib.findPi2Index(i, i, i, j, "C");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./2.));

                  col = sslib.findPi2Index(i, i, i, j, "D");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./2.));

                  // permute pop in right locus
                  col = sslib.findPi2Index(i, i, j, i, "A");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./2.));

                  col = sslib.findPi2Index(i, i, j, i, "B");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./2.));

                  col = sslib.findPi2Index(i, i, j, i, "C");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./2.));

                  col = sslib.findPi2Index(i, i, j, i, "D");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./2.));

                  // switch right appendix to j
                  col = sslib.findPi2Index(i, i, j, j, "A");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1.));

                  col = sslib.findPi2Index(i, i, j, j, "B");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1.));

                  col = sslib.findPi2Index(i, i, j, j, "C");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1.));

                  col = sslib.findPi2Index(i, i, j, j, "D");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1.));

                  // switch left appendix to j
                  col = sslib.findPi2Index(i, j, i, j, "A");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./4.));

                  col = sslib.findPi2Index(i, j, i, j, "B");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./4.));

                  col = sslib.findPi2Index(i, j, i, j, "C");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./4.));

                  col = sslib.findPi2Index(i, j, i, j, "D");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./4.));

                  // permute pop in right locus
                  col = sslib.findPi2Index(i, j, j, i, "A");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./4.));

                  col = sslib.findPi2Index(i, j, j, i, "B");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./4.));

                  col = sslib.findPi2Index(i, j, j, i, "C");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./4.));

                  col = sslib.findPi2Index(i, j, j, i, "D");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./4.));

                  // permute pops in left locus
                  col = sslib.findPi2Index(j, i, i, j, "A");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./4.));

                  col = sslib.findPi2Index(j, i, i, j, "B");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./4.));

                  col = sslib.findPi2Index(j, i, i, j, "C");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./4.));

                  col = sslib.findPi2Index(j, i, i, j, "D");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./4.));

                  // permute pop in both loci
                  col = sslib.findPi2Index(j, i, j, i, "A");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./4.));

                  col = sslib.findPi2Index(j, i, j, i, "B");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./4.));

                  col = sslib.findPi2Index(j, i, j, i, "C");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./4.));

                  col = sslib.findPi2Index(j, i, j, i, "D");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./4.));

                  // switch both apendices
                  col = sslib.findPi2Index(i, j, j, j, "A");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./2.));

                  col = sslib.findPi2Index(i, j, j, j, "B");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./2.));

                  col = sslib.findPi2Index(i, j, j, j, "C");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./2.));

                  col = sslib.findPi2Index(i, j, j, j, "D");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./2.));

                  // permute pop in left locus
                  col = sslib.findPi2Index(j, i, j, j, "A");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./2.));

                  col = sslib.findPi2Index(j, i, j, j, "B");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./2.));

                  col = sslib.findPi2Index(j, i, j, j, "C");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./2.));

                  col = sslib.findPi2Index(j, i, j, j, "D");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./2.));
                }
              }

              else if(p2 == j)
              {
                if(p3 == i) // D_j_z_ji
                {
                  // the Dz cols
                  col = sslib.findDzIndex(i, j, i);
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));

                  col = sslib.findDzIndex(j, i, i);
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));

                  // the pi2 cols
                  // imagine starting with pop indices p2 and p3 on each side of ';' character in pi2(**;**)
                  // now append i to left of p2 and right of p3

                  // as is
                  col = sslib.findPi2Index(i, j, i, i, "A");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./2.));

                  col = sslib.findPi2Index(i, j, i, i, "B");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./2.));

                  col = sslib.findPi2Index(i, j, i, i, "C");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./2.));

                  col = sslib.findPi2Index(i, j, i, i, "D");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./2.));

                  // permute pop in left locus
                  col = sslib.findPi2Index(j, i, i, i, "A");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./2.));

                  col = sslib.findPi2Index(j, i, i, i, "B");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./2.));

                  col = sslib.findPi2Index(j, i, i, i, "C");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./2.));

                  col = sslib.findPi2Index(j, i, i, i, "D");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./2.));

                  // switch right appendix to j
                  col = sslib.findPi2Index(i, j, i, j, "A");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./4.));

                  col = sslib.findPi2Index(i, j, i, j, "B");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./4.));

                  col = sslib.findPi2Index(i, j, i, j, "C");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./4.));

                  col = sslib.findPi2Index(i, j, i, j, "D");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./4.));

                  // permute pop in right locus
                  col = sslib.findPi2Index(i, j, j, i, "A");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./4.));

                  col = sslib.findPi2Index(i, j, j, i, "B");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./4.));

                  col = sslib.findPi2Index(i, j, j, i, "C");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./4.));

                  col = sslib.findPi2Index(i, j, j, i, "D");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./4.));

                  // permute pop in left locus
                  col = sslib.findPi2Index(j, i, i, j, "A");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./4.));

                  col = sslib.findPi2Index(j, i, i, j, "B");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./4.));

                  col = sslib.findPi2Index(j, i, i, j, "C");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./4.));

                  col = sslib.findPi2Index(j, i, i, j, "D");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./4.));

                  // permute pops in both loci
                  col = sslib.findPi2Index(j, i, j, i, "A");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./4.));

                  col = sslib.findPi2Index(j, i, j, i, "B");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./4.));

                  col = sslib.findPi2Index(j, i, j, i, "C");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./4.));

                  col = sslib.findPi2Index(j, i, j, i, "D");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./4.));

                  // switch left appendix to j
                  col = sslib.findPi2Index(j, j, i, i, "A");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1.));

                  col = sslib.findPi2Index(j, j, i, i, "B");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1.));

                  col = sslib.findPi2Index(j, j, i, i, "C");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1.));

                  col = sslib.findPi2Index(j, j, i, i, "D");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1.));

                  // switch both appendices
                  col = sslib.findPi2Index(j, j, i, j, "A");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./2.));

                  col = sslib.findPi2Index(j, j, i, j, "B");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./2.));

                  col = sslib.findPi2Index(j, j, i, j, "C");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./2.));

                  col = sslib.findPi2Index(j, j, i, j, "D");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./2.));

                  // permute pop in right locus
                  col = sslib.findPi2Index(j, j, j, i, "A");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./2.));

                  col = sslib.findPi2Index(j, j, j, i, "B");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./2.));

                  col = sslib.findPi2Index(j, j, j, i, "C");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./2.));

                  col = sslib.findPi2Index(j, j, j, i, "D");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./2.));
                }

                else if(p3 == j) // D_j_z_jj
                {
                  // the Dz cols
                  col = sslib.findDzIndex(i, j, j);
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));

                  col = sslib.findDzIndex(j, i, j);
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));

                  col = sslib.findDzIndex(j, j, i);
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));

                  // the pi2 cols
                  // imagine starting with pop indices p2 and p3 on each side of ';' character in pi2(**;**)
                  // now append i to left of both p2 and p3

                  // as is
                  col = sslib.findPi2Index(i, j, i, j, "A");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./4.));

                  col = sslib.findPi2Index(i, j, i, j, "B");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./4.));

                  col = sslib.findPi2Index(i, j, i, j, "C");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./4.));

                  col = sslib.findPi2Index(i, j, i, j, "D");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./4.));

                  // permute pop in left locus
                  col = sslib.findPi2Index(j, i, i, j, "A");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./4.));

                  col = sslib.findPi2Index(j, i, i, j, "B");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./4.));

                  col = sslib.findPi2Index(j, i, i, j, "C");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./4.));

                  col = sslib.findPi2Index(j, i, i, j, "D");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./4.));

                  // permute pop in right locus
                  col = sslib.findPi2Index(i, j, j, i, "A");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./4.));

                  col = sslib.findPi2Index(i, j, j, i, "B");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./4.));

                  col = sslib.findPi2Index(i, j, j, i, "C");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./4.));

                  col = sslib.findPi2Index(i, j, j, i, "D");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./4.));

                  // permute pops in both loci
                  col = sslib.findPi2Index(j, i, j, i, "A");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./4.));

                  col = sslib.findPi2Index(j, i, j, i, "B");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./4.));

                  col = sslib.findPi2Index(j, i, j, i, "C");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./4.));

                  col = sslib.findPi2Index(j, i, j, i, "D");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./4.));

                  // switch right appendix to j
                  col = sslib.findPi2Index(i, j, j, j, "A");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./2.));

                  col = sslib.findPi2Index(i, j, j, j, "B");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./2.));

                  col = sslib.findPi2Index(i, j, j, j, "C");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./2.));

                  col = sslib.findPi2Index(i, j, j, j, "D");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./2.));

                  // permute pop in left locus
                  col = sslib.findPi2Index(j, i, j, j, "A");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./2.));

                  col = sslib.findPi2Index(j, i, j, j, "B");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./2.));

                  col = sslib.findPi2Index(j, i, j, j, "C");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./2.));

                  col = sslib.findPi2Index(j, i, j, j, "D");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./2.));

                  // switch left appendix to j
                  col = sslib.findPi2Index(j, j, i, j, "A");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./2.));

                  col = sslib.findPi2Index(j, j, i, j, "B");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./2.));

                  col = sslib.findPi2Index(j, j, i, j, "C");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./2.));

                  col = sslib.findPi2Index(j, j, i, j, "D");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./2.));

                  // permute pop in right locus
                  col = sslib.findPi2Index(j, j, j, i, "A");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./2.));

                  col = sslib.findPi2Index(j, j, j, i, "B");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./2.));

                  col = sslib.findPi2Index(j, j, j, i, "C");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./2.));

                  col = sslib.findPi2Index(j, j, j, i, "D");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1./2.));

                  // switch both appendices
                  col = sslib.findPi2Index(j, j, j, j, "A");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));

                  col = sslib.findPi2Index(j, j, j, j, "B");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));

                  col = sslib.findPi2Index(j, j, j, j, "C");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));

                  col = sslib.findPi2Index(j, j, j, j, "D");
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));
                }
              }
            }
          }

          else if((*it)->getPrefix() == "pi2")
          {
            int countLeft = 0; // count of j before ';' in pi2(**;**)
            int countRight = 0; // count of j after ';' in pi2(**;**)

            if((*it)->getPopIndices()[0] == j)
              ++countLeft;

            if((*it)->getPopIndices()[1] == j)
              ++countLeft;

            if((*it)->getPopIndices()[2] == j)
              ++countRight;

            if((*it)->getPopIndices()[3] == j)
              ++countRight;

            if(countLeft + countRight > 0)  // not to populate the sparse matrix with unnecessary zeros
              coeffs.emplace_back(Eigen::Triplet<double>(row, row, -(countLeft + countRight)));

            // updates off-diagonal entries from left perspective
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

            if(countLeft > 0)  // not to populate the sparse matrix with unnecessary zeros
            {
              col = sslib.findPi2Index(p1, p2, p3, p4, "A");
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, countLeft / 16.));

              col = sslib.findPi2Index(p1, p2, p3, p4, "B");
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, countLeft / 16.));

              col = sslib.findPi2Index(p1, p2, p3, p4, "C");
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, countLeft / 16.));

              col = sslib.findPi2Index(p1, p2, p3, p4, "D");
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, countLeft / 16.));

              // permutes pops in right locus
              col = sslib.findPi2Index(p1, p2, p4, p3, "A");
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, countLeft / 16.));

              col = sslib.findPi2Index(p1, p2, p4, p3, "B");
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, countLeft / 16.));

              col = sslib.findPi2Index(p1, p2, p4, p3, "C");
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, countLeft / 16.));

              col = sslib.findPi2Index(p1, p2, p4, p3, "D");
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, countLeft / 16.));

              // permutes pops in left locus
              col = sslib.findPi2Index(p2, p1, p3, p4, "A");
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, countLeft / 16.));

              col = sslib.findPi2Index(p2, p1, p3, p4, "B");
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, countLeft / 16.));

              col = sslib.findPi2Index(p2, p1, p3, p4, "C");
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, countLeft / 16.));

              col = sslib.findPi2Index(p2, p1, p3, p4, "D");
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, countLeft / 16.));

              // permutes pops in both loci
              col = sslib.findPi2Index(p2, p1, p4, p3, "A");
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, countLeft / 16.));

              col = sslib.findPi2Index(p2, p1, p4, p3, "B");
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, countLeft / 16.));

              col = sslib.findPi2Index(p2, p1, p4, p3, "C");
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, countLeft / 16.));

              col = sslib.findPi2Index(p2, p1, p4, p3, "D");
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, countLeft / 16.));
            }

            // updates off-diagonal entries from right perspective
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

            if(countRight > 0)  // not to populate the sparse matrix with unnecessary zeros
            {
              col = sslib.findPi2Index(p1, p2, p3, p4, "A");
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, countRight / 16.));

              col = sslib.findPi2Index(p1, p2, p3, p4, "B");
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, countRight / 16.));

              col = sslib.findPi2Index(p1, p2, p3, p4, "C");
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, countRight / 16.));

              col = sslib.findPi2Index(p1, p2, p3, p4, "D");
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, countRight / 16.));

               // permutes pops in left locus
              col = sslib.findPi2Index(p2, p1, p3, p4, "A");
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, countRight / 16.));

              col = sslib.findPi2Index(p2, p1, p3, p4, "B");
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, countRight / 16.));

              col = sslib.findPi2Index(p2, p1, p3, p4, "C");
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, countRight / 16.));

              col = sslib.findPi2Index(p2, p1, p3, p4, "D");
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, countRight / 16.));

              // permutes pops in right locus
              col = sslib.findPi2Index(p1, p2, p4, p3, "A");
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, countRight / 16.));

              col = sslib.findPi2Index(p1, p2, p4, p3, "B");
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, countRight / 16.));

              col = sslib.findPi2Index(p1, p2, p4, p3, "C");
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, countRight / 16.));

              col = sslib.findPi2Index(p1, p2, p4, p3, "D");
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, countRight / 16.));

              // permutes pops in both loci
              col = sslib.findPi2Index(p2, p1, p4, p3, "A");
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, countRight / 16.));

              col = sslib.findPi2Index(p2, p1, p4, p3, "B");
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, countRight / 16.));

              col = sslib.findPi2Index(p2, p1, p4, p3, "C");
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, countRight / 16.));

              col = sslib.findPi2Index(p2, p1, p4, p3, "D");
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, countRight / 16.));
            }
          }

          else if((*it)->getPrefix() == "H")
          {
            if(childPopIdCount > 0)  // not to populate the sparse matrix with unnecessary zeros
              coeffs.emplace_back(Eigen::Triplet<double>(row, row, -childPopIdCount)); // main diagonal

            size_t p1 = (*it)->getPopIndices()[0];
            size_t p2 = (*it)->getPopIndices()[1];

            if(childPopIdCount == 1)
            {
              if(p1 == i || p2 == i) // H_ij or H_ji
              {
                col = sslib.findHetIndex(i, i, "A");
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./2.));

                col = sslib.findHetIndex(i, i, "B");
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./2.));
              }
            }

            else if(childPopIdCount == 2)
            {
              col = sslib.findHetIndex(i, j, "A");
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./2.));

              col = sslib.findHetIndex(i, j, "B");
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./2.));

              col = sslib.findHetIndex(j, i, "A");
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./2.));

              col = sslib.findHetIndex(j, i, "B");
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1./2.));
            }
          }

          else if((*it)->getPrefix() != "I")
            throw bpp::Exception("Migration::mis-specified Moment prefix: " + (*it)->getPrefix());
        }

        Eigen::SparseMatrix<double> mat(numStats, numStats);
        mat.setFromTriplets(std::begin(coeffs), std::end(coeffs));
        mat.makeCompressed();
        mat *= getParameterValue("m_" + bpp::TextTools::toString(i) + "_" + bpp::TextTools::toString(j));
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
        paramName = "m_" + bpp::TextTools::toString(i) + "_" + bpp::TextTools::toString(j);

        double prevVal = prevParams_.getParameterValue(paramName);
        double newVal = getParameterValue(paramName);

        matrices_[index] *= (newVal / prevVal);

        ++index;
      }
    }
  }

  prevParams_.matchParametersValues(getParameters());
}

