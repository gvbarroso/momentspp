/*
 * Authors: Gustavo V. Barroso
 * Created: 21/04/2023
 * Last modified: 11/05/2023
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
      std::cout << "pulse from pop " << ancFromId << " to pop " << ancToId << "; f = " << f << ", 1 - f = " << g << "\n";

      std::vector<Eigen::Triplet<double>> coeffs(0);
      coeffs.reserve(3 * sizeOfBasis);

      for(auto it = std::begin(sslib.getMoments()); it != std::end(sslib.getMoments()); ++it)
      {
        if((*it)->countInstances(ancToId) > 0) // admixture pulse is directional
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
          if((*it)->getPrefix() == "DD")
          {
            std::cout << "focal mom : " << (*it)->getName() << std::endl;

            size_t p1 = (*it)->getPopIndices()[0];
            if(p1 == ancToId)
              p1 = (*it)->getPopIndices()[1];

            // reference moments to help track which other Dz/pi2 moments contribute to focal DD
            std::shared_ptr<DzMoment> syntheticDz = std::make_shared<DzMoment>("Dz_" + bpp::TextTools::toString(p1) + "_" + bpp::TextTools::toString(ancToId) + "_" + bpp::TextTools::toString(ancToId), 0.);
            std::shared_ptr<Pi2Moment> syntheticPi2 = std::make_shared<Pi2Moment>("pi2_" + bpp::TextTools::toString(ancToId) + "_" + bpp::TextTools::toString(ancToId) + "_" + bpp::TextTools::toString(ancToId) + "_" + bpp::TextTools::toString(ancToId), 0., nullptr, nullptr);

            for(auto itCmp = std::begin(sslib.getMoments()); itCmp != std::end(sslib.getMoments()); ++itCmp)
            {
              // contributions from Dz moments
              if(syntheticDz->isAdmixAdjacent(*itCmp, ancFromId, ancToId))
              {
                col = sslib.findCompressedIndex(*itCmp);

                size_t t = 1 + (*it)->countInstances(ancToId);
                size_t a = 1 + ((*itCmp)->getPopIndices()[0] == ancToId);
                size_t b = t - a;
                double x = (0.25 * (*it)->countInstances(ancToId) * std::pow(-1., 4. - (*itCmp)->countInstances(ancFromId)) * std::pow(g, a) * std::pow(f, b)) / ((*it)->getNumberOfAliases() + 1);

                coeffs.emplace_back(Eigen::Triplet<double>(row, col, x));
              }

              // contributions from pi2 moments
              else if(syntheticPi2->isAdmixAdjacent(*itCmp, ancFromId, ancToId) && (*it)->countInstances(ancToId) == 2)
              {
                col = sslib.findCompressedIndex(*itCmp);

                size_t parentPopIdCount =(*itCmp)->countInstances(ancFromId);
                double x = std::pow(g, 2.) * std::pow(f, 2.);
                double y = std::pow(-1., 4. - parentPopIdCount) * x;
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, y));
              }
            }
          }

          else if((*it)->getPrefix() == "Dz")
          {
            if((*it)->getPopIndices()[0] == ancToId)
            {
              if((*it)->countInstances(ancToId) == 3)
              {
                //std::cout << "focal Dz: " << (*it)->getName() << std::endl;

                // a reference pi2 moment to help track which other pi2 moments contribute to focal Dz
                std::shared_ptr<Pi2Moment> synthetic = std::make_shared<Pi2Moment>("pi2_" + bpp::TextTools::toString(ancToId) + "_" + bpp::TextTools::toString(ancToId) + "_" + bpp::TextTools::toString(ancToId) + "_" + bpp::TextTools::toString(ancToId), 0., nullptr, nullptr);

                std::cout << "synthetic: " << synthetic->getName() << std::endl;
                for(auto itPi2 = std::begin(sslib.getMoments()); itPi2 != std::end(sslib.getMoments()); ++itPi2)
                {
                  if(synthetic->isAdmixAdjacent((*itPi2), ancFromId, ancToId))
                  {
                    auto tmp = std::dynamic_pointer_cast<Pi2Moment>(*itPi2);
                    std::cout << "adjacent: " << tmp->getName() << std::endl;

                    col = sslib.findCompressedIndex(tmp);

                    size_t parentPopIdCount = tmp->countInstances(ancFromId);
                    size_t childPopIdCount = tmp->countInstances(ancToId);

                    size_t l = tmp->getLeftHetStat()->isCrossPop();
                    size_t r = tmp->getRightHetStat()->isCrossPop();
                    size_t z = l + r;
                    size_t w1 = parentPopIdCount > childPopIdCount;
                    size_t w2 = childPopIdCount > parentPopIdCount;
                    size_t x = 1 + tmp->hasPopIndex(ancFromId) + w1 - z;
                    size_t y = 1 + tmp->hasPopIndex(ancToId) + w2 - z;

                    double c = std::pow(f, x) * std::pow(1. - f, y) * std::pow(1. - 2. * f, z);
                    double d = std::pow(-1., y - 1);

                    coeffs.emplace_back(Eigen::Triplet<double>(row, col, d * c));
                  }
                }
              }

              else {
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
