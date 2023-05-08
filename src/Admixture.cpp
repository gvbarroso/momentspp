/*
 * Authors: Gustavo V. Barroso
 * Created: 21/04/2023
 * Last modified: 08/05/2023
 *
 */

#include <typeinfo>

#include "Admixture.hpp"

void Admixture::setUpMatrices_(const SumStatsLibrary& sslib)
{
  size_t numPops = popIndices_.size();
  size_t sizeOfBasis = sslib.getSizeOfBasis();
  matrices_.reserve(getParameters().size()); // max. one matrix per population

  for(size_t i = 0; i < numPops; ++i)
  {
    size_t ancFromId = 0;
    size_t ancToId = 0;
    double f = -1.;
    double g = -1.;

    for(size_t j = 0; j < numPops; ++j)
    {
      if(littleAdmixMat_(i, j) > 0.)
      {
        ancToId = popIndices_[i];
        ancFromId = popIndices_[j];
        f = littleAdmixMat_(i, j);
        g = 1. - f;
      }
    }

    if(f > 0.)
    {
      //std::cout << "pulse from pop " << ancFromId << " to pop " << ancToId << "; f = " << f << ", 1 - f = " << g << "\n";

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
            if((*it)->isAdmixAdjacent(*it2nd, ancFromId, ancToId))
            {
              col = sslib.findCompressedIndex(*it2nd);

              double mig = (*it)->countInstances(ancToId) - (*it2nd)->countInstances(ancToId);
              double nat = (*it2nd)->countInstances(ancToId);
              double y = std::pow(g, nat) * std::pow(f, mig) / ((*it)->getNumberOfAliases() + 1);

              coeffs.emplace_back(Eigen::Triplet<double>(row, col, y));
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

            double c = (f * g * (static_cast<double>(p2 == p3) / 2. - 0.25)) / ((*it)->getNumberOfAliases() + 1);
            col = sslib.findCompressedIndex(sslib.findDzIndex(p1, p2, p3));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, c));
            col = sslib.findCompressedIndex(sslib.findDzIndex(p1, p3, p2));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, c));

            p2 = ancToId;

            c = f * g * (static_cast<double>(p2 == p3) / 2. - 0.25) / ((*it)->getNumberOfAliases() + 1);
            col = sslib.findCompressedIndex(sslib.findDzIndex(p1, p2, p3));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, c));

            p2 = ancFromId;
            p3 = ancFromId;

            c = f * g * (static_cast<double>(p2 == p3) / 2. - 0.25) / ((*it)->getNumberOfAliases() + 1);
            col = sslib.findCompressedIndex(sslib.findDzIndex(p1, p2, p3));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, c));

            if((*it)->countInstances(ancToId) == 2)
            {
              p1 = ancFromId;
              p2 = ancFromId;
              p3 = ancToId;

              c = (f * g * (static_cast<double>(p2 == p3) / 2. - 0.25)) / ((*it)->getNumberOfAliases() + 1);
              col = sslib.findCompressedIndex(sslib.findDzIndex(p1, p2, p3));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, c));
              col = sslib.findCompressedIndex(sslib.findDzIndex(p1, p3, p2));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, c));

              p2 = ancToId;

              c = f * g * (static_cast<double>(p2 == p3) / 2. - 0.25) / ((*it)->getNumberOfAliases() + 1);
              col = sslib.findCompressedIndex(sslib.findDzIndex(p1, p2, p3));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, c));

              p2 = ancFromId;
              p3 = ancFromId;

              c = f * g * (static_cast<double>(p2 == p3) / 2. - 0.25) / ((*it)->getNumberOfAliases() + 1);
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
                  double x = std::pow(g, childPopIdCount) * std::pow(f, parentPopIdCount);
                  double y = std::pow(-1., 4. - parentPopIdCount) * x;
                  coeffs.emplace_back(Eigen::Triplet<double>(row, col, y));
                }
              }
            }
          }

          else if((*it)->getPrefix() == "Dz")
          {
            if((*it)->getPopIndices()[0] == ancToId)
            {
              std::vector<size_t> popIds = { (*it)->getPopIndices()[1], ancFromId, (*it)->getPopIndices()[2], ancFromId };

              size_t refCount = std::count(std::begin(popIds), std::end(popIds), ancToId);
              size_t count = std::count(std::begin(popIds), std::end(popIds), ancToId);
              double y = f * g * std::pow(-1., count - refCount);

              col = sslib.findCompressedIndex(sslib.findPi2Index(popIds[0], popIds[1], popIds[2], popIds[3]));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, y));
              col = sslib.findCompressedIndex(sslib.findPi2Index(popIds[0], popIds[1], popIds[3], popIds[2]));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, y));
              col = sslib.findCompressedIndex(sslib.findPi2Index(popIds[1], popIds[0], popIds[2], popIds[3]));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, y));
              col = sslib.findCompressedIndex(sslib.findPi2Index(popIds[1], popIds[0], popIds[3], popIds[2]));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, y));

              popIds[3] = ancToId; // switch right appendix to childPopId
              count = std::count(std::begin(popIds), std::end(popIds), ancToId);
              y = f * g * std::pow(-1., count - refCount);

              col = sslib.findCompressedIndex(sslib.findPi2Index(popIds[0], popIds[1], popIds[2], popIds[3]));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, y));
              col = sslib.findCompressedIndex(sslib.findPi2Index(popIds[0], popIds[1], popIds[3], popIds[2]));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, y));
              col = sslib.findCompressedIndex(sslib.findPi2Index(popIds[1], popIds[0], popIds[2], popIds[3]));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, y));
              col = sslib.findCompressedIndex(sslib.findPi2Index(popIds[1], popIds[0], popIds[3], popIds[2]));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, y));

              popIds[3] = ancFromId; // back
              popIds[1] = ancToId; // switch left appendix to childPopId
              count = std::count(std::begin(popIds), std::end(popIds), ancToId);
              y = f * g * std::pow(-1., count - refCount);

              col = sslib.findCompressedIndex(sslib.findPi2Index(popIds[0], popIds[1], popIds[2], popIds[3]));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, y));
              col = sslib.findCompressedIndex(sslib.findPi2Index(popIds[0], popIds[1], popIds[3], popIds[2]));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, y));
              col = sslib.findCompressedIndex(sslib.findPi2Index(popIds[1], popIds[0], popIds[2], popIds[3]));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, y));
              col = sslib.findCompressedIndex(sslib.findPi2Index(popIds[1], popIds[0], popIds[3], popIds[2]));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, y));

              popIds[3] = ancToId; // have both switched
              count = std::count(std::begin(popIds), std::end(popIds), ancToId);
              y = f * g * std::pow(-1., count - refCount);

              col = sslib.findCompressedIndex(sslib.findPi2Index(popIds[0], popIds[1], popIds[2], popIds[3]));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, y));
              col = sslib.findCompressedIndex(sslib.findPi2Index(popIds[0], popIds[1], popIds[3], popIds[2]));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, y));
              col = sslib.findCompressedIndex(sslib.findPi2Index(popIds[1], popIds[0], popIds[2], popIds[3]));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, y));
              col = sslib.findCompressedIndex(sslib.findPi2Index(popIds[1], popIds[0], popIds[3], popIds[2]));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, y));
            }
          }
        }
      } // ends loop over moments

      Eigen::SparseMatrix<double> mat(sizeOfBasis, sizeOfBasis);
      mat.setFromTriplets(std::begin(coeffs), std::end(coeffs));
      mat.makeCompressed();
      matrices_.emplace_back(mat);
    }
  }

  setIdentity_(sizeOfBasis);
  assembleTransitionMatrix_();
}

void Admixture::updateMatrices_()
{
  size_t numPops = littleAdmixMat_.innerSize();
  std::string paramName = "";

  for(size_t i = 0; i < numPops; ++i)
  {
    size_t id = popIndices_[i];

    for(size_t j = 0; j < numPops; ++j)
    {
      size_t jd = popIndices_[j];

      if(id != jd)
      {
        if(littleAdmixMat_(i, j) > 0.)
        {
          paramName = "a_" + bpp::TextTools::toString(id) + "_" + bpp::TextTools::toString(jd);

          if(hasParameter(paramName))
          {
            double newVal = getParameterValue(paramName);
            littleAdmixMat_(i, j) = newVal;
          }

          else
          {
            paramName = "a_" + bpp::TextTools::toString(jd) + "_" + bpp::TextTools::toString(id); // switch order
            double newVal = getParameterValue(paramName);
            littleAdmixMat_(i, j) = 1. - newVal;
          }
        }
      }
    }
  }

  //setUpMatrices_(sslib); // WARNING
  assembleTransitionMatrix_();
  prevParams_.matchParametersValues(getParameters());
}
