/*
 * Authors: Gustavo V. Barroso
 * Created: 21/04/2023
 * Last modified: 02/05/2023
 *
 */

#include <typeinfo>

#include "Admixture.hpp"

void Admixture::setUpMatrices_(const SumStatsLibrary& sslib)
{
  size_t numPops = popIndices_.size();
  size_t sizeOfBasis = sslib.getSizeOfBasis();
  matrices_.reserve(getParameters().size());

  for(size_t i = 0; i < numPops; ++i)
  {
    size_t ancFromId = popIndices_[i]; // contributes f

    for(size_t j = 0; j < numPops; ++j)
    {
      size_t ancToId = popIndices_[j]; // contributes 1-f
      double f = littleAdmixMat_(i, j);

      if(ancFromId != ancToId && f > 0.)
      {
        std::vector<Eigen::Triplet<double>> coeffs(0);
        coeffs.reserve(3 * sizeOfBasis);

        for(auto it = std::begin(sslib.getMoments()); it != std::end(sslib.getMoments()); ++it)
        {
          if((*it)->countInstances(ancToId) > 0) // admixture is directional
          {
            int row = sslib.findCompressedIndex(*it);
            int col = -1; // inits column index to out-of-bounds

            // contributions from moments of the same prefix
            for(auto it2nd = std::begin(sslib.getMoments()); it2nd != std::end(sslib.getMoments()); ++it2nd)
            {
              if((*it)->getPrefix() == (*it2nd)->getPrefix())
              {
                col = sslib.findCompressedIndex(*it2nd);

                size_t parentPopIdCount =(*it2nd)->countInstances(ancFromId);
                size_t childPopIdCount = (*it2nd)->countInstances(ancToId);
                double y = std::pow(1. - f, childPopIdCount) * std::pow(f, parentPopIdCount);
                double z = ((*it)->isAdmixAdjacent(*it2nd, ancFromId, ancToId) * y) / ((*it)->getNumberOfAliases() + 1);

                coeffs.emplace_back(Eigen::Triplet<double>(row, col, z));
              }
            }

            // contributions from moments of distinct prefixes
            if((*it)->getPrefix() == "DD" && (*it)->hasPopIndex(ancToId))
            {
              size_t p1 = (*it)->getPopIndices()[0];
              if(p1 == ancToId)
                p1 = (*it)->getPopIndices()[1];

              size_t p2 = ancFromId;
              size_t p3 = ancToId;

              double c = (f * (1. - f) * (static_cast<double>(p2 == p3) / 2. - 0.25)) / ((*it)->getNumberOfAliases() + 1);
              col = sslib.findCompressedIndex(sslib.findDzIndex(p1, p2, p3));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, c));
              col = sslib.findCompressedIndex(sslib.findDzIndex(p1, p3, p2));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, c));

              p2 = ancToId;

              c = f * (1. - f) * (static_cast<double>(p2 == p3) / 2. - 0.25) / ((*it)->getNumberOfAliases() + 1);
              col = sslib.findCompressedIndex(sslib.findDzIndex(p1, p2, p3));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, c));

              p2 = ancFromId;
              p3 = ancFromId;

              c = f * (1. - f) * (static_cast<double>(p2 == p3) / 2. - 0.25) / ((*it)->getNumberOfAliases() + 1);
              col = sslib.findCompressedIndex(sslib.findDzIndex(p1, p2, p3));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, c));

              if((*it)->countInstances(ancToId) == 2)
              {
                p1 = ancFromId;
                p2 = ancFromId;
                p3 = ancToId;

                c = (f * (1. - f) * (static_cast<double>(p2 == p3) / 2. - 0.25)) / ((*it)->getNumberOfAliases() + 1);
                col = sslib.findCompressedIndex(sslib.findDzIndex(p1, p2, p3));
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, c));
                col = sslib.findCompressedIndex(sslib.findDzIndex(p1, p3, p2));
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, c));

                p2 = ancToId;

                c = f * (1. - f) * (static_cast<double>(p2 == p3) / 2. - 0.25) / ((*it)->getNumberOfAliases() + 1);
                col = sslib.findCompressedIndex(sslib.findDzIndex(p1, p2, p3));
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, c));

                p2 = ancFromId;
                p3 = ancFromId;

                c = f * (1. - f) * (static_cast<double>(p2 == p3) / 2. - 0.25) / ((*it)->getNumberOfAliases() + 1);
                col = sslib.findCompressedIndex(sslib.findDzIndex(p1, p2, p3));
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, c));

                // contributions from pi2 moments
                for(auto itPi2 = std::begin(sslib.getMoments()); itPi2 != std::end(sslib.getMoments()); ++itPi2)
                {
                  if((*itPi2)->getPrefix() == "pi2")
                  {
                    col = sslib.findCompressedIndex(*itPi2);

                    size_t parentPopIdCount =(*itPi2)->countInstances(ancFromId);
                    size_t childPopIdCount = (*itPi2)->countInstances(ancToId);
                    double x = std::pow(1. - f, childPopIdCount) * std::pow(f, parentPopIdCount);
                    double y = std::pow(-1., 4. - parentPopIdCount) * x;
                    coeffs.emplace_back(Eigen::Triplet<double>(row, col, y));
                  }
                }
              }
            }
          /*
          else if((*it)->getPrefix() == "Dz")
          {
            std::vector<size_t> popIds = (*it)->getPopIndices();

            for(size_t l = 0; l < popIds.size(); ++ l) // contributions from the Dz cols
            {
              if(popIds[l] == jd) // if entry matches childPopId
              {
                popIds[l] = id; // assign to focal parentPopId
                col = sslib.findCompressedIndex(sslib.findDzIndex(popIds[0], popIds[1], popIds[2]));
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, 0.5));
                col = sslib.findCompressedIndex(sslib.findDzIndex(popIds[0], popIds[2], popIds[1]));
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, 0.5));
                popIds[l] = jd; // recycle
              }
            }

            if((*it)->getPopIndices()[0] == jd) // contributions from pi2 moments
            {
              // imagine starting with pop indices p2 and p3 on each side of ';' character in pi2(**;**)
              // append parentPopId (i) to the right of both p2 and p3
              // find pi2 statistics by left-right permuting + replacing appendixes

              popIds.clear();
              popIds.resize(4);

              popIds[0] = (*it)->getPopIndices()[1];
              popIds[1] = id;
              popIds[2] = (*it)->getPopIndices()[2];
              popIds[3] = id;

              size_t refCount = std::count(std::begin(popIds), std::end(popIds), jd);
              size_t count = std::count(std::begin(popIds), std::end(popIds), jd);
              double f = std::pow(-1., count - refCount);

              col = sslib.findCompressedIndex(sslib.findPi2Index(popIds[0], popIds[1], popIds[2], popIds[3]));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, f));
              col = sslib.findCompressedIndex(sslib.findPi2Index(popIds[0], popIds[1], popIds[3], popIds[2]));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, f));
              col = sslib.findCompressedIndex(sslib.findPi2Index(popIds[1], popIds[0], popIds[2], popIds[3]));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, f));
              col = sslib.findCompressedIndex(sslib.findPi2Index(popIds[1], popIds[0], popIds[3], popIds[2]));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, f));

              popIds[3] = jd; // switch right appendix to childPopId
              count = std::count(std::begin(popIds), std::end(popIds), jd);
              f = std::pow(-1., count - refCount);

              col = sslib.findCompressedIndex(sslib.findPi2Index(popIds[0], popIds[1], popIds[2], popIds[3]));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, f));
              col = sslib.findCompressedIndex(sslib.findPi2Index(popIds[0], popIds[1], popIds[3], popIds[2]));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, f));
              col = sslib.findCompressedIndex(sslib.findPi2Index(popIds[1], popIds[0], popIds[2], popIds[3]));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, f));
              col = sslib.findCompressedIndex(sslib.findPi2Index(popIds[1], popIds[0], popIds[3], popIds[2]));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, f));

              popIds[3] = id; // back
              popIds[1] = jd; // switch left appendix to childPopId
              count = std::count(std::begin(popIds), std::end(popIds), jd);
              f = std::pow(-1., count - refCount);

              col = sslib.findCompressedIndex(sslib.findPi2Index(popIds[0], popIds[1], popIds[2], popIds[3]));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, f));
              col = sslib.findCompressedIndex(sslib.findPi2Index(popIds[0], popIds[1], popIds[3], popIds[2]));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, f));
              col = sslib.findCompressedIndex(sslib.findPi2Index(popIds[1], popIds[0], popIds[2], popIds[3]));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, f));
              col = sslib.findCompressedIndex(sslib.findPi2Index(popIds[1], popIds[0], popIds[3], popIds[2]));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, f));

              popIds[3] = jd; // have both switched
              count = std::count(std::begin(popIds), std::end(popIds), jd);
              f = std::pow(-1., count - refCount);

              col = sslib.findCompressedIndex(sslib.findPi2Index(popIds[0], popIds[1], popIds[2], popIds[3]));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, f));
              col = sslib.findCompressedIndex(sslib.findPi2Index(popIds[0], popIds[1], popIds[3], popIds[2]));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, f));
              col = sslib.findCompressedIndex(sslib.findPi2Index(popIds[1], popIds[0], popIds[2], popIds[3]));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, f));
              col = sslib.findCompressedIndex(sslib.findPi2Index(popIds[1], popIds[0], popIds[3], popIds[2]));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, f));
            }
          }

          else if((*it)->getPrefix() != "I")
            throw bpp::Exception("Admixture::mis-specified Moment prefix: " + (*it)->getPrefix());
          */
          }
        }

        Eigen::SparseMatrix<double> mat(sizeOfBasis, sizeOfBasis);
        mat.setFromTriplets(std::begin(coeffs), std::end(coeffs));
        mat.makeCompressed();
        matrices_.emplace_back(mat);
      }
    }
  }

  setIdentity_(sizeOfBasis);
  assembleTransitionMatrix_();
}

void Admixture::updateMatrices_()
{
  size_t numPops = littleAdmixMat_.innerSize();
  size_t index = 0;
  std::string paramName = "";

  for(size_t i = 0; i < numPops; ++i)
  {
    size_t id = popIndices_[i];

    for(size_t j = 0; j < numPops; ++j)
    {
      size_t jd = popIndices_[j];

      if(id != jd)
      {
        paramName = "m_" + bpp::TextTools::toString(id) + "_" + bpp::TextTools::toString(jd);

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
