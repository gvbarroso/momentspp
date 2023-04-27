/*
 * Authors: Gustavo V. Barroso
 * Created: 21/04/2023
 * Last modified: 27/04/2023
 *
 */

#include <typeinfo>

#include "Admixture.hpp"

void Admixture::setUpMatrices_(const SumStatsLibrary& sslib)
{
  size_t numPops = popIndices_.size();
  size_t sizeOfBasis = sslib.getSizeOfBasis();
  matrices_.reserve(numPops * (numPops - 1));

  for(size_t i = 0; i < numPops; ++i)
  {
    size_t ancOneId = popIndices_[i]; // contributes f

    for(size_t j = 0; j < numPops; ++j)
    {
      size_t ancTwoId = popIndices_[j]; // contributes 1-f

      if(ancOneId != ancTwoId)
      {
        double f = littleAdmixMat_(i, j);

        std::vector<Eigen::Triplet<double>> coeffs(0);
        coeffs.reserve(3 * sizeOfBasis);

        for(auto it = std::begin(sslib.getBasis()); it != std::end(sslib.getBasis()); ++it)
        {
          int row = it - std::begin(sslib.getBasis()); // row index
          int col = -1; // inits column index to out-of-bounds

          int childPopIdCount = static_cast<int>((*it)->countInstances(ancTwoId));
          int parentPopIdCount = static_cast<int>((*it)->countInstances(ancOneId));

          double x = std::pow(f, childPopIdCount);
          double y = std::pow(1. - f, parentPopIdCount);

          // constributions from moments of the same kind
          coeffs.emplace_back(Eigen::Triplet<double>(row, row, x * y));

          for(auto it2nd = std::begin(sslib.getBasis()); it2nd != std::end(sslib.getBasis()); ++it2nd)
          {
            if(typeid((*it).get()) == typeid((*it2nd).get()))
            {
              col = it2nd - std::begin(sslib.getBasis());
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, adjacencyMat_(row, col) * x * y));
            }
          }

          /*
          if((*it)->getPrefix() == "DD")
          {
            std::vector<size_t> popIds = (*it)->getPopIndices();

            for(size_t l = 0; l < popIds.size(); ++ l) // contributions from the DD cols
            {
              if(popIds[l] == jd) // if entry matches childPopId
              {
                popIds[l] = id; // assign to focal parentPopId
                col = sslib.findCompressedIndex(sslib.findDdIndex(popIds[0], popIds[1]));
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, 0.5));
                col = sslib.findCompressedIndex(sslib.findDdIndex(popIds[1], popIds[0]));
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, 0.5));
                popIds[l] = jd; // recycle
              }
            }

            if((*it)->hasPopIndex(j))
            {
              size_t p1 = popIds[0];
              if(p1 == jd)
                p1 = popIds[1];

              size_t p2 = id;
              size_t p3 = jd;

              double f = static_cast<double>(p2 == p3) / 2. - 0.25;
              col = sslib.findCompressedIndex(sslib.findDzIndex(p1, p2, p3));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, childPopIdCount * f));
              col = sslib.findCompressedIndex(sslib.findDzIndex(p1, p3, p2));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, childPopIdCount * f));

              p2 = jd;

              f = static_cast<double>(p2 == p3) / 2. - 0.25;
              col = sslib.findCompressedIndex(sslib.findDzIndex(p1, p2, p3));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, childPopIdCount * f));

              p2 = id;
              p3 = id;

              f = static_cast<double>(p2 == p3) / 2. - 0.25;
              col = sslib.findCompressedIndex(sslib.findDzIndex(p1, p2, p3));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, childPopIdCount * f));
            }
          }

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

        Eigen::SparseMatrix<double> mat(sizeOfBasis, sizeOfBasis);
        mat.setFromTriplets(std::begin(coeffs), std::end(coeffs));
        mat.makeCompressed();
        std::cout << mat << std::endl;
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

void Admixture::setUpAdjacencyMatrix_(const SumStatsLibrary& sslib)
{
  adjacencyMat_.resize(sslib.getSizeOfBasis(), sslib.getSizeOfBasis());
  adjacencyMat_.setZero();

  for(size_t i = 0; i < sslib.getSizeOfBasis(); ++i)
  {
    for(size_t j = 0; j < sslib.getSizeOfBasis(); ++j)
    {
      if(sslib.getBasis()[i]->getPrefix() == sslib.getBasis()[j]->getPrefix())
      {
        if(sslib.getBasis()[i]->isAdjacent(sslib.getBasis()[j]))
          adjacencyMat_(i, j) = 1.;
      }
    }
  }
}
