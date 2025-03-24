/*
 * Authors: Gustavo V. Barroso
 * Created: 21/03/2025
 * Last modified: 25/03/2025
 *
 */

#include <typeinfo>

#include "NeutralAdmixture.hpp"

void NeutralAdmixture::setUpMatrices_(const SumStatsLibrary& sslib)
{
  size_t numPops = popIndices_.size();
  size_t sizeOfBasis = sslib.getSizeOfBasis();
  matrices_.reserve(getParameters().size()); // max. one matrix per population

  for(size_t i = 0; i < numPops; ++i)
  {
    size_t ancFromId = 0;
    size_t ancToId = 0;
    long double f = -1.;

    for(size_t j = 0; j < numPops; ++j)
    {
      if(littleAdmixMat_(i, j) > 0.)
      {
        ancFromId = popIndices_[i];
        ancToId = popIndices_[j];
        f = littleAdmixMat_(i, j);
      }
    }

    if(f > 0.)
    {
      std::vector<Eigen::Triplet<long double>> coeffs(0);
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

              long double mig = (*it)->countInstances(ancToId) - (*it2nd)->countInstances(ancToId);
              long double nat = (*it2nd)->countInstances(ancToId);
              long double y = std::pow((1. - f), nat) * std::pow(f, mig) / ((*it)->getNumberOfAliases() + 1);

              coeffs.emplace_back(Eigen::Triplet<long double>(row, col, y));
            }
          }

          // contributions from moments of distinct prefixes
          if((*it)->getPrefix() == "DD")
          {
            size_t p1 = (*it)->getPopIndices()[0];

            if(p1 == ancToId)
              p1 = (*it)->getPopIndices()[1];

            // reference moments to help track which other Dr/pi2 moments contribute to focal DD
            std::shared_ptr<DrMoment> syntheticDr = std::make_shared<DrMoment>("Dr_" + bpp::TextTools::toString(p1) + "_" + bpp::TextTools::toString(ancToId) + "_" + bpp::TextTools::toString(ancToId), 0.);
            std::shared_ptr<Pi2Moment> syntheticPi2 = std::make_shared<Pi2Moment>("pi2_" + bpp::TextTools::toString(ancToId) + "_" + bpp::TextTools::toString(ancToId) + "_" + bpp::TextTools::toString(ancToId) + "_" + bpp::TextTools::toString(ancToId), 0., nullptr, nullptr);

            for(auto itCmp = std::begin(sslib.getMoments()); itCmp != std::end(sslib.getMoments()); ++itCmp)
            {
              // contributions from Dr moments
              if(syntheticDr->isAdmixAdjacent(*itCmp, ancFromId, ancToId))
              {
                col = sslib.findCompressedIndex(*itCmp);

                size_t t = 1 + (*it)->countInstances(ancToId);
                size_t a = 1 + ((*itCmp)->getPopIndices()[0] == ancToId);
                size_t b = t - a;
                long double c = std::pow(-1., (*itCmp)->countInstances(ancFromId));

                if((*itCmp)->getPopIndices()[0] == ancFromId)
                  c *= -1.;

                long double x = (0.25 * (*it)->countInstances(ancToId) * c * std::pow((1. - f), a) * std::pow(f, b)) / ((*it)->getNumberOfAliases() + 1);

                coeffs.emplace_back(Eigen::Triplet<long double>(row, col, x));
              }

              // contributions from pi2 moments
              else if(syntheticPi2->isAdmixAdjacent(*itCmp, ancFromId, ancToId) && (*it)->countInstances(ancToId) == 2)
              {
                col = sslib.findCompressedIndex(*itCmp);
                size_t parentPopIdCount =(*itCmp)->countInstances(ancFromId);
                long double x = std::pow((1. - f), 2.) * std::pow(f, 2.);
                long double y = std::pow(-1., parentPopIdCount) * x;
                coeffs.emplace_back(Eigen::Triplet<long double>(row, col, y));
              }
            }
          }

          else if((*it)->getPrefix() == "Dr")
          {
            if((*it)->getPopIndices()[0] == ancToId) // has pi2 contributions
            {
              // a reference pi2 moment to help track which other pi2 moments contribute to focal Dr
              std::shared_ptr<Pi2Moment> syntheticPi2 = std::make_shared<Pi2Moment>("pi2_" + bpp::TextTools::toString(ancToId) + "_" + bpp::TextTools::toString((*it)->getPopIndices()[1]) + "_" + bpp::TextTools::toString(ancToId) + "_" + bpp::TextTools::toString((*it)->getPopIndices()[2]), 0., nullptr, nullptr);

              // case 1: Dr TTT
              if((*it)->countInstances(ancToId) == 3)
              {
                col = sslib.findCompressedIndex(sslib.findPi2Index(ancFromId, ancFromId, ancFromId, ancFromId));
                coeffs.emplace_back(Eigen::Triplet<long double>(row, col, 4. * std::pow(f, 3.) * (1. - f)));

                col = sslib.findCompressedIndex(sslib.findPi2Index(ancFromId, ancToId, ancFromId, ancFromId));
                coeffs.emplace_back(Eigen::Triplet<long double>(row, col, 2. * std::pow(f, 2.) * (1. - f) * (1. - 2 * f)));
                col = sslib.findCompressedIndex(sslib.findPi2Index(ancToId, ancFromId, ancFromId, ancFromId));
                coeffs.emplace_back(Eigen::Triplet<long double>(row, col, 2. * std::pow(f, 2.) * (1. - f) * (1. - 2 * f)));
                col = sslib.findCompressedIndex(sslib.findPi2Index(ancFromId, ancFromId, ancToId, ancFromId));
                coeffs.emplace_back(Eigen::Triplet<long double>(row, col, 2. * std::pow(f, 2.) * (1. - f) * (1. - 2 * f)));
                col = sslib.findCompressedIndex(sslib.findPi2Index(ancFromId, ancFromId, ancFromId, ancToId));
                coeffs.emplace_back(Eigen::Triplet<long double>(row, col, 2. * std::pow(f, 2.) * (1. - f) * (1. - 2 * f)));

                col = sslib.findCompressedIndex(sslib.findPi2Index(ancFromId, ancFromId, ancToId, ancToId));
                coeffs.emplace_back(Eigen::Triplet<long double>(row, col, -4. * std::pow(f, 2.) * std::pow(1. - f, 2.)));
                col = sslib.findCompressedIndex(sslib.findPi2Index(ancToId, ancToId, ancFromId, ancFromId));
                coeffs.emplace_back(Eigen::Triplet<long double>(row, col, -4. * std::pow(f, 2.) * std::pow(1. - f, 2.)));

                col = sslib.findCompressedIndex(sslib.findPi2Index(ancFromId, ancToId, ancFromId, ancToId));
                coeffs.emplace_back(Eigen::Triplet<long double>(row, col, f * (1. - f) * std::pow(1. - 2. * f, 2.)));
                col = sslib.findCompressedIndex(sslib.findPi2Index(ancToId, ancFromId, ancFromId, ancToId));
                coeffs.emplace_back(Eigen::Triplet<long double>(row, col, f * (1. - f) * std::pow(1. - 2. * f, 2.)));
                col = sslib.findCompressedIndex(sslib.findPi2Index(ancFromId, ancToId, ancToId, ancFromId));
                coeffs.emplace_back(Eigen::Triplet<long double>(row, col, f * (1. - f) * std::pow(1. - 2. * f, 2.)));
                col = sslib.findCompressedIndex(sslib.findPi2Index(ancToId, ancFromId, ancToId, ancFromId));
                coeffs.emplace_back(Eigen::Triplet<long double>(row, col, f * (1. - f) * std::pow(1. - 2. * f, 2.)));

                col = sslib.findCompressedIndex(sslib.findPi2Index(ancToId, ancToId, ancFromId, ancToId));
                coeffs.emplace_back(Eigen::Triplet<long double>(row, col, -2. * f * std::pow(1. - f, 2.) * (1. - 2. * f)));
                col = sslib.findCompressedIndex(sslib.findPi2Index(ancToId, ancToId, ancToId, ancFromId));
                coeffs.emplace_back(Eigen::Triplet<long double>(row, col, -2. * f * std::pow(1. - f, 2.) * (1. - 2. * f)));
                col = sslib.findCompressedIndex(sslib.findPi2Index(ancFromId, ancToId, ancToId, ancToId));
                coeffs.emplace_back(Eigen::Triplet<long double>(row, col, -2. * f * std::pow(1. - f, 2.) * (1. - 2. * f)));
                col = sslib.findCompressedIndex(sslib.findPi2Index(ancToId, ancFromId, ancToId, ancToId));
                coeffs.emplace_back(Eigen::Triplet<long double>(row, col, -2. * f * std::pow(1. - f, 2.) * (1. - 2. * f)));

                col = sslib.findCompressedIndex(sslib.findPi2Index(ancToId, ancToId, ancToId, ancToId));
                coeffs.emplace_back(Eigen::Triplet<long double>(row, col, 4. * f * std::pow(1. - f, 3.)));
              }

              // case 2: Dr TFT / TTF NOTE "forced" through some differences with the "naive" Mathematica notebook to match moments.LD admix matrix
              else if(((*it)->getPopIndices()[1] == ancToId && (*it)->getPopIndices()[2] == ancFromId) || ((*it)->getPopIndices()[1] == ancFromId && (*it)->getPopIndices()[2] == ancToId))
              {
                col = sslib.findCompressedIndex(sslib.findPi2Index(ancFromId, ancFromId, ancFromId, ancFromId));
                coeffs.emplace_back(Eigen::Triplet<long double>(row, col, 2. * std::pow(f, 2.) * (1. - f)));

                col = sslib.findCompressedIndex(sslib.findPi2Index(ancFromId, ancToId, ancFromId, ancFromId));
                coeffs.emplace_back(Eigen::Triplet<long double>(row, col, 2. * f * (1. - f) * (1. - 3. * f)));

                col = sslib.findCompressedIndex(sslib.findPi2Index(ancFromId, ancFromId, ancToId, ancToId));
                coeffs.emplace_back(Eigen::Triplet<long double>(row, col, -2. * f * std::pow(1. - f, 2.)));

                col = sslib.findCompressedIndex(sslib.findPi2Index(ancFromId, ancToId, ancFromId, ancToId));
                coeffs.emplace_back(Eigen::Triplet<long double>(row, col, -2. * f * (1. - f) * (1. - 2. * f)));

                col = sslib.findCompressedIndex(sslib.findPi2Index(ancToId, ancToId, ancFromId, ancToId));
                coeffs.emplace_back(Eigen::Triplet<long double>(row, col, 2. * f * std::pow(1. - f, 2.)));
              }

              // case 3: Dr TFF
              else if((*it)->getPopIndices()[1] == ancFromId && (*it)->getPopIndices()[2] == ancFromId)
              {
                col = sslib.findCompressedIndex(sslib.findPi2Index(ancFromId, ancFromId, ancFromId, ancFromId));
                coeffs.emplace_back(Eigen::Triplet<long double>(row, col, 4. * f * (1. - f)));

                col = sslib.findCompressedIndex(sslib.findPi2Index(ancFromId, ancToId, ancFromId, ancFromId));
                coeffs.emplace_back(Eigen::Triplet<long double>(row, col, -2. * f * (1. - f)));
                col = sslib.findCompressedIndex(sslib.findPi2Index(ancToId, ancFromId, ancFromId, ancFromId));
                coeffs.emplace_back(Eigen::Triplet<long double>(row, col, -2. * f * (1. - f)));
                col = sslib.findCompressedIndex(sslib.findPi2Index(ancFromId, ancFromId, ancToId, ancFromId));
                coeffs.emplace_back(Eigen::Triplet<long double>(row, col, -2. * f * (1. - f)));
                col = sslib.findCompressedIndex(sslib.findPi2Index(ancFromId, ancFromId, ancFromId, ancToId));
                coeffs.emplace_back(Eigen::Triplet<long double>(row, col, -2. * f * (1. - f)));

                col = sslib.findCompressedIndex(sslib.findPi2Index(ancFromId, ancToId, ancFromId, ancToId));
                coeffs.emplace_back(Eigen::Triplet<long double>(row, col, f * (1. - f)));
                col = sslib.findCompressedIndex(sslib.findPi2Index(ancFromId, ancToId, ancToId, ancFromId));
                coeffs.emplace_back(Eigen::Triplet<long double>(row, col, f * (1. - f)));
                col = sslib.findCompressedIndex(sslib.findPi2Index(ancToId, ancFromId, ancFromId, ancToId));
                coeffs.emplace_back(Eigen::Triplet<long double>(row, col, f * (1. - f)));
                col = sslib.findCompressedIndex(sslib.findPi2Index(ancToId, ancFromId, ancToId, ancFromId));
                coeffs.emplace_back(Eigen::Triplet<long double>(row, col, f * (1. - f)));
              }

              // case 4: Dr TTx / TxT
              else if((*it)->getPopIndices()[1] == ancToId || (*it)->getPopIndices()[2] == ancToId)
              {
                for(auto itCmp = std::begin(sslib.getMoments()); itCmp != std::end(sslib.getMoments()); ++itCmp)
                {
                  if(syntheticPi2->isAdmixAdjacent(*itCmp, ancFromId, ancToId))
                  {
                    auto tmp = std::dynamic_pointer_cast<Pi2Moment>(*itCmp);
                    col = sslib.findCompressedIndex(tmp);

                    long double c = std::pow(-1., tmp->countInstances(ancToId) + 1);
                    long double z = (tmp->getLeftHetStat()->isCrossPop() && tmp->getRightHetStat()->isCrossPop());
                    long double x = 1 + (tmp->getLeftHetStat()->countInstances(ancFromId) == 2 || tmp->getRightHetStat()->countInstances(ancFromId) == 2);
                    long double y = 3. - z - x;
                    long double w = std::pow(f, x) * std::pow((1. - f), y) * std::pow(1. - 2 * f, z);

                    if(z == 0)
                    {
                      w *= 2.;
                      if(x == 2)
                        w *= -1;
                    }

                    coeffs.emplace_back(Eigen::Triplet<long double>(row, col, c * w));
                  }
                }
              }

              // case 5: Dr TFx / TxF
              else if((*it)->getPopIndices()[1] == ancFromId || (*it)->getPopIndices()[2] == ancFromId)
              {
                for(auto itCmp = std::begin(sslib.getMoments()); itCmp != std::end(sslib.getMoments()); ++itCmp)
                {
                  if(syntheticPi2->isAdmixAdjacent(*itCmp, ancFromId, ancToId))
                  {
                    auto tmp = std::dynamic_pointer_cast<Pi2Moment>(*itCmp);
                    col = sslib.findCompressedIndex(tmp);
                    long double c = std::pow(-1., tmp->countInstances(ancToId));
                    coeffs.emplace_back(Eigen::Triplet<long double>(row, col, c * 2 * f * (1. - f)));
                  }
                }
              }

              // case 6: Dr Txx / Txx
              else
              {
                for(auto itCmp = std::begin(sslib.getMoments()); itCmp != std::end(sslib.getMoments()); ++itCmp)
                {
                  if(syntheticPi2->isAdmixAdjacent(*itCmp, ancFromId, ancToId))
                  {
                    auto tmp = std::dynamic_pointer_cast<Pi2Moment>(*itCmp);
                    col = sslib.findCompressedIndex(tmp);
                    long double c = std::pow(-1., tmp->countInstances(ancToId));
                    long double d = 2. * (1. + ((*it)->getPopIndices()[1] == (*it)->getPopIndices()[2]));
                    coeffs.emplace_back(Eigen::Triplet<long double>(row, col, c * d * f * (1. - f)));
                  }
                }
              }
            }
          }
        }
      } // ends loop over moments

      Eigen::SparseMatrix<long double> mat(sizeOfBasis, sizeOfBasis);
      mat.setFromTriplets(std::begin(coeffs), std::end(coeffs));
      mat.makeCompressed();
      matrices_.emplace_back(mat);
    }
  }

  assembleTransitionMatrix_();
}

void NeutralAdmixture::updateMatrices_()
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
            long double newVal = getParameterValue(paramName);
            littleAdmixMat_(i, j) = newVal;
          }

          else
          {
            paramName = "a_" + bpp::TextTools::toString(jd) + "_" + bpp::TextTools::toString(id); // switch order
            long double newVal = getParameterValue(paramName);
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

// overrides AbstractOperator because NeutralAdmixture works differently
void NeutralAdmixture::assembleTransitionMatrix_()
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
    long double rowSum = 0.;

    for(int j = 0; j < transition_.cols(); ++j)
      rowSum += transition_.coeffRef(i, j);

    if(rowSum == 0.)
      transition_.coeffRef(i, i) = 1.;
  }
}
