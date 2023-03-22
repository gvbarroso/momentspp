/*
 * Authors: Gustavo V. Barroso
 * Created: 10/08/2022
 * Last modified: 21/03/2023
 *
 */


#include "Migration.hpp"


void Migration::setUpMatrices_(const SumStatsLibrary& sslib)
{
  // m_ij is the forward migration rate from pop i to pop j (backwards, the prob that lineage in j has parent in i)
  size_t numPops = fetchNumPops();
  size_t sizeOfBasis = sslib.getSizeOfBasis();
  matrices_.reserve(numPops * (numPops - 1));

  for(size_t i = 0; i < numPops; ++i) // for i in m_ij (i => parent pop ID)
  {
    for(size_t j = 0; j < numPops; ++j) // for j in m_ij (j => child pop ID)
    {
      if(i != j)
      {
        std::vector<Eigen::Triplet<double>> coeffs(0);
        coeffs.reserve(3 * sizeOfBasis);

        // although we have pi2(i,j;k,l) moments, we need only to loop over pops twice because we take care of the k,l pop indices in this loop over stats
        for(auto it = std::begin(sslib.getBasis()); it != std::end(sslib.getBasis()); ++it)
        {
          int row = it - std::begin(sslib.getBasis()); // row index
          int col = -1; // inits column index to out-of-bounds
          int childPopIdCount = static_cast<int>((*it)->countInstances(j));

          if(childPopIdCount != 0) // not to populate the sparse matrix with unnecessary zeros
            coeffs.emplace_back(Eigen::Triplet<double>(row, row, -childPopIdCount));

          if((*it)->getPrefix() == "DD")
          {
            if(childPopIdCount == 1)
            {
              col = sslib.findCompressedIndex(sslib.findDdIndex(i, i));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));
            }

            else if(childPopIdCount == 2)
            {
              col = sslib.findCompressedIndex(sslib.findDdIndex(i, j));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));

              col = sslib.findCompressedIndex(sslib.findDdIndex(j, i));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));
            }

            for(size_t k = sslib.getNumDDStats(); k < (sslib.getNumDDStats() + sslib.getNumDzStats()); ++k)
            {
              assert(sslib.getMoments()[k]->getPrefix() == "Dz");

              size_t pop1 = sslib.getMoments()[k]->getPopIndices()[0];
              size_t pop2 = sslib.getMoments()[k]->getPopIndices()[1];
              size_t pop3 = sslib.getMoments()[k]->getPopIndices()[2];

              if(((pop1 == i && childPopIdCount == 1) || (pop1 == j && childPopIdCount == 2)) && ((pop2 == i || pop2 == j) && (pop3 == i || pop3 == j)))
              {
                double f = static_cast<double>(pop2 == pop3) / 2. - 0.25;
                col = sslib.findCompressedIndex(k);
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, childPopIdCount * f));
              }
            }
          }

          else if((*it)->getPrefix() == "Dz")
          {
            size_t p1 = (*it)->getPopIndices()[0]; // D pop
            size_t p2 = (*it)->getPopIndices()[1]; // first z pop
            size_t p3 = (*it)->getPopIndices()[2]; // second z pop

            if(p1 == i) // D_i_z_**
            {
              if((p2 == i && p3 == j) || (p2 == j && p3 == i)) // D_i_z_ij || D_i_z_ji
              {
                col = sslib.findCompressedIndex(sslib.findDzIndex(i, i, i));
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));
              }

              else if(p2 == j && p3 == j) // D_i_z_jj
              {
                col = sslib.findCompressedIndex(sslib.findDzIndex(i, i, j));
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));

                col = sslib.findCompressedIndex(sslib.findDzIndex(i, j, i));
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
                  col = sslib.findCompressedIndex(sslib.findDzIndex(i, i, i));
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));

                  // the pi2 cols (minding permutations of derived/ancestral encoded in the suffixes, see Pi2Moment class)
                  // imagine starting with pop indices p2 and p3 on each side of ';' character in pi2(**;**) -> [i*;*i]
                  // now append i to the right of both p2 and p3

                  // as is
                  col = sslib.findCompressedIndex(sslib.findPi2Index(i, i, i, i));
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 4.));

                  // switch right appendix to j
                  col = sslib.findCompressedIndex(sslib.findPi2Index(i, i, i, j));
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -2.));

                  // permute pop in right locus
                  col = sslib.findCompressedIndex(sslib.findPi2Index(i, i, j, i));
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -2.));

                  // switch left appendix to j
                  col = sslib.findCompressedIndex(sslib.findPi2Index(i, j, i, i));
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -2.));

                  // permute pop in left locus
                  col = sslib.findCompressedIndex(sslib.findPi2Index(j, i, i, i));
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -2.));

                  // switch both
                  col = sslib.findCompressedIndex(sslib.findPi2Index(i, j, i, j));
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));

                  // permute pop in right locus
                  col = sslib.findCompressedIndex(sslib.findPi2Index(i, j, j, i));
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));

                  // permute pop in left locus
                  col = sslib.findCompressedIndex(sslib.findPi2Index(j, i, i, j));
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));

                  // permute pops in both loci
                  col = sslib.findCompressedIndex(sslib.findPi2Index(j, i, j, i));
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));
                }

                else if(p3 == j) // D_j_z_ij
                {
                  // the Dz cols
                  col = sslib.findCompressedIndex(sslib.findDzIndex(j, i, i));
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));

                  col = sslib.findCompressedIndex(sslib.findDzIndex(i, i, j));
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 0.5));

                  col = sslib.findCompressedIndex(sslib.findDzIndex(i, j, i));
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 0.5));

                  // the pi2 cols
                  // imagine starting with pop indices p2 and p3 on each side of ';' character in pi2(**;**) -> [i*;*j]
                  // now append i to right the of p2 and to the left of p3

                  // as is
                  col = sslib.findCompressedIndex(sslib.findPi2Index(i, i, i, j));
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 2.));

                  // permute pop in right locus
                  col = sslib.findCompressedIndex(sslib.findPi2Index(i, i, j, i));
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 2.));

                  // switch right appendix to j
                  col = sslib.findCompressedIndex(sslib.findPi2Index(i, i, j, j));
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -4.));

                  // switch left appendix to j
                  col = sslib.findCompressedIndex(sslib.findPi2Index(i, j, i, j));
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1.));

                  // permute pop in right locus
                  col = sslib.findCompressedIndex(sslib.findPi2Index(i, j, j, i));
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1.));

                  // permute pops in left locus
                  col = sslib.findCompressedIndex(sslib.findPi2Index(j, i, i, j));
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1.));

                  // permute pop in both loci
                  col = sslib.findCompressedIndex(sslib.findPi2Index(j, i, j, i));
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1.));

                  // switch both apendices
                  col = sslib.findCompressedIndex(sslib.findPi2Index(i, j, j, j));
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 2.));

                  // permute pop in left locus
                  col = sslib.findCompressedIndex(sslib.findPi2Index(j, i, j, j));
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 2.));
                }
              }

              else if(p2 == j)
              {
                if(p3 == i) // D_j_z_ji
                {
                  // the Dz cols
                  col = sslib.findCompressedIndex(sslib.findDzIndex(j, i, i));
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));

                  col = sslib.findCompressedIndex(sslib.findDzIndex(i, j, i));
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 0.5));

                  col = sslib.findCompressedIndex(sslib.findDzIndex(i, i, j));
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 0.5));

                  // the pi2 cols
                  // imagine starting with pop indices p2 and p3 on each side of ';' character in pi2(**;**)
                  // now append i to left of p2 and right of p3

                  // as is
                  col = sslib.findCompressedIndex(sslib.findPi2Index(i, j, i, i));
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 2.));

                  // permute pop in left locus
                  col = sslib.findCompressedIndex(sslib.findPi2Index(j, i, i, i));
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 2.));

                  // switch right appendix to j
                  col = sslib.findCompressedIndex(sslib.findPi2Index(i, j, i, j));
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1.));

                  // permute pop in right locus
                  col = sslib.findCompressedIndex(sslib.findPi2Index(i, j, j, i));
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1.));

                  // permute pop in left locus
                  col = sslib.findCompressedIndex(sslib.findPi2Index(j, i, i, j));
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1.));

                  // permute pops in both loci
                  col = sslib.findCompressedIndex(sslib.findPi2Index(j, i, j, i));
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -1.));

                  // switch left appendix to j
                  col = sslib.findCompressedIndex(sslib.findPi2Index(j, j, i, i));
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -4.));

                  // switch both appendices
                  col = sslib.findCompressedIndex(sslib.findPi2Index(j, j, i, j));
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 2.));

                  // permute pop in right locus
                  col = sslib.findCompressedIndex(sslib.findPi2Index(j, j, j, i));
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 2.));
                }

                else if(p3 == j) // D_j_z_jj
                {
                  // the Dz cols
                  col = sslib.findCompressedIndex(sslib.findDzIndex(i, j, j));
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));

                  col = sslib.findCompressedIndex(sslib.findDzIndex(j, i, j));
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));

                  col = sslib.findCompressedIndex(sslib.findDzIndex(j, j, i));
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));

                  // the pi2 cols
                  // imagine starting with pop indices p2 and p3 on each side of ';' character in pi2(**;**)
                  // now append i to left of both p2 and p3

                  // as is
                  col = sslib.findCompressedIndex(sslib.findPi2Index(i, j, i, j));
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));

                  // permute pop in left locus
                  col = sslib.findCompressedIndex(sslib.findPi2Index(j, i, i, j));
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));

                  // permute pop in right locus
                  col = sslib.findCompressedIndex(sslib.findPi2Index(i, j, j, i));
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));

                  // permute pops in both loci
                  col = sslib.findCompressedIndex(sslib.findPi2Index(j, i, j, i));
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));

                  // switch right appendix to j
                  col = sslib.findCompressedIndex(sslib.findPi2Index(i, j, j, j));
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -2.));

                  // permute pop in left locus
                  col = sslib.findCompressedIndex(sslib.findPi2Index(j, i, j, j));
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -2.));

                  // switch left appendix to j
                  col = sslib.findCompressedIndex(sslib.findPi2Index(j, j, i, j));
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -2.));

                  // permute pop in right locus
                  col = sslib.findCompressedIndex(sslib.findPi2Index(j, j, j, i));
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, -2.));

                  // switch both appendices
                  col = sslib.findCompressedIndex(sslib.findPi2Index(j, j, j, j));
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, 4.));
                }
              }
            }
          }

          else if((*it)->getPrefix() == "pi2") // finds pi2 moments that are a single migration away (1-hop neighbors)
          {
            std::vector<size_t> popIds = (*it)->getPopIndices();

            for(size_t l = 0; l < popIds.size(); ++ l) // l -> 0:3
            {
              if(popIds[l] == j) // if entry matches childPopId
              {
                popIds[l] = i; // assign to focal parentPopId
                col = sslib.findCompressedIndex(sslib.findPi2Index(popIds[0], popIds[1], popIds[2], popIds[3]));
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));
                popIds[l] = j; // recycle
              }
            }
          }

          else if((*it)->getPrefix() == "H")
          {
            std::vector<size_t> popIds = (*it)->getPopIndices();

            for(size_t l = 0; l < popIds.size(); ++ l) // l -> 0:1
            {
              if(popIds[l] == j) // if entry matches childPopId
              {
                popIds[l] = i; // assign to focal parentPopId
                col = sslib.findCompressedIndex(sslib.findHetIndex(popIds[0], popIds[1]));
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));
                popIds[l] = j; // recycle
              }
            }
          }

          else if((*it)->getPrefix() != "I")
            throw bpp::Exception("Migration::mis-specified Moment prefix: " + (*it)->getPrefix());
        }

        Eigen::SparseMatrix<double> mat(sizeOfBasis, sizeOfBasis);
        mat.setFromTriplets(std::begin(coeffs), std::end(coeffs));
        mat.makeCompressed();
        mat *= getParameterValue("m_" + bpp::TextTools::toString(i) + "_" + bpp::TextTools::toString(j));
        matrices_.emplace_back(mat);
      } // ends if(i != j)
    }
  }

  setIdentity_(sizeOfBasis);
  assembleTransitionMatrix_();
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

        if(newVal != prevVal)
          matrices_[index] *= (newVal / prevVal);

        ++index;
      }
    }
  }

  assembleTransitionMatrix_();
  prevParams_.matchParametersValues(getParameters());
}

