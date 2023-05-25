/*
 * Authors: Gustavo V. Barroso
 * Created: 21/04/2023
 * Last modified: 25/05/2023
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
    size_t ancTwoId = 0;
    size_t ancOneId = 0;
    double f = -1.;

    for(size_t j = 0; j < numPops; ++j)
    {
      if(littleAdmixMat_(i, j) > 0.)
      {
        ancOneId = popIndices_[i];
        ancTwoId = popIndices_[j];
        f = littleAdmixMat_(i, j);

        //std::cout << "from " << ancOneId << " to " << ancTwoId << ", f = " << f << std::endl;
        //std::cout << littleAdmixMat_ << std::endl;
      }
    }

    if(f > 0.)
    {
      std::vector<Eigen::Triplet<double>> coeffs(0);
      coeffs.reserve(3 * sizeOfBasis);

      for(auto it = std::begin(sslib.getMoments()); it != std::end(sslib.getMoments()); ++it)
      {
        if((*it)->countInstances(ancOneId) > 0)
        {
          int row = sslib.findCompressedIndex(*it);
          int col = -1; // inits column index to out-of-bounds

          // contributions from moments of the same prefix
          for(auto it2nd = std::begin(sslib.getMoments()); it2nd != std::end(sslib.getMoments()); ++it2nd)
          {
            if((*it)->isAdmixAdjacent(*it2nd, ancTwoId, ancOneId))
            {
              col = sslib.findCompressedIndex(*it2nd);

              double mig = (*it)->countInstances(ancOneId) - (*it2nd)->countInstances(ancOneId);
              double nat = (*it2nd)->countInstances(ancOneId);
              double y = std::pow((1. - f), nat) * std::pow(f, mig) / ((*it)->getNumberOfAliases() + 1);

              coeffs.emplace_back(Eigen::Triplet<double>(row, col, y));
            }
          }

          // contributions from moments of distinct prefixes
          if((*it)->getPrefix() == "DD")
          {
            size_t p1 = (*it)->getPopIndices()[0];

            if(p1 == ancOneId)
              p1 = (*it)->getPopIndices()[1];

            // reference moments to help track which other Dz/pi2 moments contribute to focal DD
            std::shared_ptr<DzMoment> syntheticDz = std::make_shared<DzMoment>("Dz_" + bpp::TextTools::toString(p1) + "_" + bpp::TextTools::toString(ancOneId) + "_" + bpp::TextTools::toString(ancOneId), 0.);
            std::shared_ptr<Pi2Moment> syntheticPi2 = std::make_shared<Pi2Moment>("pi2_" + bpp::TextTools::toString(ancOneId) + "_" + bpp::TextTools::toString(ancOneId) + "_" + bpp::TextTools::toString(ancOneId) + "_" + bpp::TextTools::toString(ancOneId), 0., nullptr, nullptr);

            for(auto itCmp = std::begin(sslib.getMoments()); itCmp != std::end(sslib.getMoments()); ++itCmp)
            {
              // contributions from Dz moments
              if(syntheticDz->isAdmixAdjacent(*itCmp, ancTwoId, ancOneId))
              {
                col = sslib.findCompressedIndex(*itCmp);

                size_t t = 1 + (*it)->countInstances(ancOneId);
                size_t a = 1 + ((*itCmp)->getPopIndices()[0] == ancOneId);
                size_t b = t - a;
                double c = std::pow(-1., (*itCmp)->countInstances(ancTwoId));

                if((*itCmp)->getPopIndices()[0] == ancTwoId)
                  c *= -1.;

                double x = (0.25 * (*it)->countInstances(ancOneId) * c * std::pow((1. - f), a) * std::pow(f, b)) / ((*it)->getNumberOfAliases() + 1);

                coeffs.emplace_back(Eigen::Triplet<double>(row, col, x));
              }

              // contributions from pi2 moments
              else if(syntheticPi2->isAdmixAdjacent(*itCmp, ancTwoId, ancOneId) && (*it)->countInstances(ancOneId) == 2)
              {
                col = sslib.findCompressedIndex(*itCmp);
                size_t parentPopIdCount =(*itCmp)->countInstances(ancTwoId);
                double x = std::pow((1. - f), 2.) * std::pow(f, 2.);
                double y = std::pow(-1., parentPopIdCount) * x;
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, y));
              }
            }
          }

          else if((*it)->getPrefix() == "Dz")
          {
            if((*it)->getPopIndices()[0] == ancOneId) // has pi2 contributions
            {
              // a reference pi2 moment to help track which other pi2 moments contribute to focal Dz
              std::shared_ptr<Pi2Moment> syntheticPi2 = std::make_shared<Pi2Moment>("pi2_" + bpp::TextTools::toString(ancOneId) + "_" + bpp::TextTools::toString((*it)->getPopIndices()[1]) + "_" + bpp::TextTools::toString(ancOneId) + "_" + bpp::TextTools::toString((*it)->getPopIndices()[2]), 0., nullptr, nullptr);

              // case 1: Dz TTT
              if((*it)->countInstances(ancOneId) == 3)
              {
                col = sslib.findCompressedIndex(sslib.findPi2Index(ancTwoId, ancTwoId, ancTwoId, ancTwoId));
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, 4. * std::pow(f, 3.) * (1. - f)));

                col = sslib.findCompressedIndex(sslib.findPi2Index(ancTwoId, ancOneId, ancTwoId, ancTwoId));
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, 2. * std::pow(f, 2.) * (1. - f) * (1. - 2 * f)));
                col = sslib.findCompressedIndex(sslib.findPi2Index(ancOneId, ancTwoId, ancTwoId, ancTwoId));
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, 2. * std::pow(f, 2.) * (1. - f) * (1. - 2 * f)));
                col = sslib.findCompressedIndex(sslib.findPi2Index(ancTwoId, ancTwoId, ancOneId, ancTwoId));
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, 2. * std::pow(f, 2.) * (1. - f) * (1. - 2 * f)));
                col = sslib.findCompressedIndex(sslib.findPi2Index(ancTwoId, ancTwoId, ancTwoId, ancOneId));
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, 2. * std::pow(f, 2.) * (1. - f) * (1. - 2 * f)));

                col = sslib.findCompressedIndex(sslib.findPi2Index(ancTwoId, ancTwoId, ancOneId, ancOneId));
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, -4. * std::pow(f, 2.) * std::pow(1. - f, 2.)));
                col = sslib.findCompressedIndex(sslib.findPi2Index(ancOneId, ancOneId, ancTwoId, ancTwoId));
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, -4. * std::pow(f, 2.) * std::pow(1. - f, 2.)));

                col = sslib.findCompressedIndex(sslib.findPi2Index(ancTwoId, ancOneId, ancTwoId, ancOneId));
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, f * (1. - f) * std::pow(1. - 2. * f, 2.)));
                col = sslib.findCompressedIndex(sslib.findPi2Index(ancOneId, ancTwoId, ancTwoId, ancOneId));
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, f * (1. - f) * std::pow(1. - 2. * f, 2.)));
                col = sslib.findCompressedIndex(sslib.findPi2Index(ancTwoId, ancOneId, ancOneId, ancTwoId));
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, f * (1. - f) * std::pow(1. - 2. * f, 2.)));
                col = sslib.findCompressedIndex(sslib.findPi2Index(ancOneId, ancTwoId, ancOneId, ancTwoId));
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, f * (1. - f) * std::pow(1. - 2. * f, 2.)));

                col = sslib.findCompressedIndex(sslib.findPi2Index(ancOneId, ancOneId, ancTwoId, ancOneId));
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, -2. * f * std::pow(1. - f, 2.) * (1. - 2. * f)));
                col = sslib.findCompressedIndex(sslib.findPi2Index(ancOneId, ancOneId, ancOneId, ancTwoId));
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, -2. * f * std::pow(1. - f, 2.) * (1. - 2. * f)));
                col = sslib.findCompressedIndex(sslib.findPi2Index(ancTwoId, ancOneId, ancOneId, ancOneId));
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, -2. * f * std::pow(1. - f, 2.) * (1. - 2. * f)));
                col = sslib.findCompressedIndex(sslib.findPi2Index(ancOneId, ancTwoId, ancOneId, ancOneId));
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, -2. * f * std::pow(1. - f, 2.) * (1. - 2. * f)));

                col = sslib.findCompressedIndex(sslib.findPi2Index(ancOneId, ancOneId, ancOneId, ancOneId));
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, 4. * f * std::pow(1. - f, 3.)));
              }
              // case 2: Dz TFT / TTF NOTE "forced" through some differences with the "naive" Mathematica notebook to match moments.LD admix matrix
              else if(((*it)->getPopIndices()[1] == ancOneId && (*it)->getPopIndices()[2] == ancTwoId) || ((*it)->getPopIndices()[1] == ancTwoId && (*it)->getPopIndices()[2] == ancOneId))
              {
                col = sslib.findCompressedIndex(sslib.findPi2Index(ancTwoId, ancTwoId, ancTwoId, ancTwoId));
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, 2. * std::pow(f, 2.) * (1. - f)));

                col = sslib.findCompressedIndex(sslib.findPi2Index(ancTwoId, ancOneId, ancTwoId, ancTwoId));
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, 2. * f * (1. - f) * (1. - 3. * f)));

                col = sslib.findCompressedIndex(sslib.findPi2Index(ancTwoId, ancTwoId, ancOneId, ancOneId));
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, -2. * f * std::pow(1. - f, 2.)));

                col = sslib.findCompressedIndex(sslib.findPi2Index(ancTwoId, ancOneId, ancTwoId, ancOneId));
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, -2. * f * (1. - f) * (1. - 2. * f)));

                col = sslib.findCompressedIndex(sslib.findPi2Index(ancOneId, ancOneId, ancTwoId, ancOneId));
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, 2. * f * std::pow(1. - f, 2.)));
              }
              // case 3: Dz TFF
              else if((*it)->getPopIndices()[1] == ancTwoId && (*it)->getPopIndices()[2] == ancTwoId)
              {
                col = sslib.findCompressedIndex(sslib.findPi2Index(ancTwoId, ancTwoId, ancTwoId, ancTwoId));
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, 4. * f * (1. - f)));

                col = sslib.findCompressedIndex(sslib.findPi2Index(ancTwoId, ancOneId, ancTwoId, ancTwoId));
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, -2. * f * (1. - f)));
                col = sslib.findCompressedIndex(sslib.findPi2Index(ancOneId, ancTwoId, ancTwoId, ancTwoId));
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, -2. * f * (1. - f)));
                col = sslib.findCompressedIndex(sslib.findPi2Index(ancTwoId, ancTwoId, ancOneId, ancTwoId));
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, -2. * f * (1. - f)));
                col = sslib.findCompressedIndex(sslib.findPi2Index(ancTwoId, ancTwoId, ancTwoId, ancOneId));
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, -2. * f * (1. - f)));

                col = sslib.findCompressedIndex(sslib.findPi2Index(ancTwoId, ancOneId, ancTwoId, ancOneId));
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, f * (1. - f)));
                col = sslib.findCompressedIndex(sslib.findPi2Index(ancTwoId, ancOneId, ancOneId, ancTwoId));
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, f * (1. - f)));
                col = sslib.findCompressedIndex(sslib.findPi2Index(ancOneId, ancTwoId, ancTwoId, ancOneId));
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, f * (1. - f)));
                col = sslib.findCompressedIndex(sslib.findPi2Index(ancOneId, ancTwoId, ancOneId, ancTwoId));
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, f * (1. - f)));
              }

              // case 4: Dz TTx / TxT
              else if((*it)->getPopIndices()[1] == ancOneId || (*it)->getPopIndices()[2] == ancOneId)
              {
                for(auto itCmp = std::begin(sslib.getMoments()); itCmp != std::end(sslib.getMoments()); ++itCmp)
                {
                  if(syntheticPi2->isAdmixAdjacent(*itCmp, ancTwoId, ancOneId))
                  {
                    auto tmp = std::dynamic_pointer_cast<Pi2Moment>(*itCmp);
                    col = sslib.findCompressedIndex(tmp);

                    double c = std::pow(-1., tmp->countInstances(ancOneId) + 1);
                    double z = (tmp->getLeftHetStat()->isCrossPop() && tmp->getRightHetStat()->isCrossPop());
                    double x = 1 + (tmp->getLeftHetStat()->countInstances(ancTwoId) == 2 || tmp->getRightHetStat()->countInstances(ancTwoId) == 2);
                    double y = 3. - z - x;
                    double w = std::pow(f, x) * std::pow((1. - f), y) * std::pow(1. - 2 * f, z);

                    if(z == 0)
                    {
                      w *= 2.;
                      if(x == 2)
                        w *= -1;
                    }

                    coeffs.emplace_back(Eigen::Triplet<double>(row, col, c * w));
                  }
                }
              }
              // case 5: Dz TFx / TxF
              else if((*it)->getPopIndices()[1] == ancTwoId || (*it)->getPopIndices()[2] == ancTwoId)
              {
                for(auto itCmp = std::begin(sslib.getMoments()); itCmp != std::end(sslib.getMoments()); ++itCmp)
                {
                  if(syntheticPi2->isAdmixAdjacent(*itCmp, ancTwoId, ancOneId))
                  {
                    auto tmp = std::dynamic_pointer_cast<Pi2Moment>(*itCmp);
                    col = sslib.findCompressedIndex(tmp);
                    double c = std::pow(-1., tmp->countInstances(ancOneId));
                    coeffs.emplace_back(Eigen::Triplet<double>(row, col, c * 2 * f * (1. - f)));
                  }
                }
              }
              // case 6: Dz Txx / Txx
              else
              {
                for(auto itCmp = std::begin(sslib.getMoments()); itCmp != std::end(sslib.getMoments()); ++itCmp)
                {
                  if(syntheticPi2->isAdmixAdjacent(*itCmp, ancTwoId, ancOneId))
                  {
                    auto tmp = std::dynamic_pointer_cast<Pi2Moment>(*itCmp);
                    col = sslib.findCompressedIndex(tmp);
                    double c = std::pow(-1., tmp->countInstances(ancOneId));
                    coeffs.emplace_back(Eigen::Triplet<double>(row, col, c * 4 * f * (1. - f)));
                  }
                }
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

  //setUpMatrices_(sslib); // WARNING inefficient
  assembleTransitionMatrix_();
  prevParams_.matchParametersValues(getParameters());
}

// overrides AbstractOperator because Admixture works differently
void Admixture::assembleTransitionMatrix_()
{
  transition_ = matrices_[0]; // inits to "delta" matrix

  if(matrices_.size() > 1)
  {
    for(size_t i = 1; i < matrices_.size(); ++i)
      transition_ += matrices_[i];
  }

  // only adds 1 to main diagonal of empty rows
  for(int i = 0; i < transition_.rows(); ++i)
  {
    double rowSum = 0.;

    for(int j = 0; j < transition_.cols(); ++j)
      rowSum += transition_.coeffRef(i, j);

    if(rowSum == 0.)
      transition_.coeffRef(i, i) = 1.;
  }
}
