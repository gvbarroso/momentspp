/*
 * Authors: Gustavo V. Barroso
 * Created: 09/08/2022
 * Last modified: 28/09/2023
 *
 */

#include <Bpp/Text/TextTools.h>

#include "Drift.hpp"

int Drift::computeDMainDiagContribution_(std::shared_ptr<Moment> mom, size_t id)
{
  size_t totalPopIdCount = mom->countInstances(id) + mom->getPopFactorPower(id);

  int a = 0;
  int b = -1;

  for(size_t i = 0; i < totalPopIdCount; ++i)
  {
    a += b;
    --b;
  }

  return a;
}

int Drift::computeDOffDiagContribution_(std::shared_ptr<Moment> mom, size_t id)
{
  size_t totalPopIdCount = mom->countInstances(id) + mom->getPopFactorPower(id);

  int a = 0;
  int b = 1;

  for(size_t i = 0; i < totalPopIdCount; ++i)
  {
    a += b;
    ++b;
  }

  return a;
}

int Drift::computeDrMainDiagContribution_(std::shared_ptr<Moment> mom, size_t id)
{
  size_t popIdCount = mom->countInstances(id);
  size_t popIdPower = mom->getPopFactorPower(id);
  size_t totalPopIdCount = popIdCount + popIdPower;

  if(popIdCount == 2)
  {
    int a = -2;
    int b = -1;

    for(size_t i = 0; i < totalPopIdCount - 1; ++i)
    {
      a += b;
      --b;
    }

    return a;
  }

  else if(popIdCount == 1)
  {
    int a = 0;
    int b = -1;

    for(size_t i = 0; i < totalPopIdCount - 1; ++i)
    {
      a += b;
      --b;
    }

    return a;
  }

  else // popIdCount == 0
    return 0;
}

int Drift::computeDDMainDiagContribution_(std::shared_ptr<Moment> mom, size_t id)
{
  size_t popIdCount = mom->countInstances(id);
  size_t popIdPower = mom->getPopFactorPower(id);
  size_t totalPopIdCount = popIdCount + popIdPower;

  if(popIdCount == 2)
  {
    int a = 0;
    int b = -3;

    for(size_t i = 0; i < totalPopIdCount - 1; ++i)
    {
      a += b;
      --b;
    }

    return a;
  }

  else if(popIdCount == 1)
  {
    int a = -2;
    int b = -1;

    for(size_t i = 0; i < totalPopIdCount - 1; ++i)
    {
      a += b;
      --b;
    }

    return a;
  }

  else // popIdCount == 0
    return 0;
}

int Drift::computePi2MainDiagContribution_(std::shared_ptr<Moment> mom, size_t id)
{
  size_t totalPopIdCount = mom->countInstances(id) + mom->getPopFactorPower(id);

  int a = -1;
  int b = -1;

  for(size_t i = 0; i < totalPopIdCount; ++i)
  {
    a += b;
    --b;
  }

  return a;
}

void Drift::setUpMatrices_(const SumStatsLibrary& sslib)
{
  size_t numPops = getParameters().size();
  size_t sizeOfBasis = sslib.getSizeOfBasis();
  matrices_.reserve(numPops);

  for(size_t i = 0; i < numPops; ++i)
  {
    size_t id = popIndices_[i];
    std::vector<Eigen::Triplet<double>> coeffs(0);
    coeffs.reserve(sizeOfBasis);

    for(auto it = std::begin(sslib.getBasis()); it != std::end(sslib.getBasis()); ++it)
    {
      int row = it - std::begin(sslib.getBasis());
      int col = -1;

      size_t popIdCount = (*it)->countInstances(id);
      int popIdPower = (*it)->getPopFactorPower(id);

      if((*it)->getPrefix() == "D")
      {
        if(popIdCount == 1)
        {
          int x = computeDMainDiagContribution_(*it, id);
          coeffs.emplace_back(Eigen::Triplet<double>(row, row, x));

          if(popIdPower > 1)
          {
            std::vector<size_t> factorIds = (*it)->getFactorIndices();
            sslib.dropFactorIds(factorIds, id, 2);

            int y = computeDOffDiagContribution_(*it, id);
            col = sslib.findCompressedIndex(sslib.getMoment("D", { id }, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, y));
          }
        }
      }

      else if((*it)->getPrefix() == "Dr")
      {
        if((*it)->getPopIndices()[0] == id) // D_id_r_*
        {
          int x = computeDrMainDiagContribution_(*it, id);
          coeffs.emplace_back(Eigen::Triplet<double>(row, row, x));

          if(popIdCount == 2) // D_id_r_id
          {
            if(popIdPower > 0)
            {
              std::vector<size_t> factorIds = (*it)->getFactorIndices();
              sslib.dropFactorIds(factorIds, id, 1);

              if((popIdPower > (sslib.getFactorOrder() + 1)) || (((*it)->getFactorPower() > sslib.getFactorOrder()) && popIdPower > 1))
                sslib.dropFactorIds(factorIds, id, 1);

              while(factorIds.size() > static_cast<size_t>(sslib.getFactorOrder())) // NOTE
                factorIds.pop_back();

              col = sslib.findCompressedIndex(sslib.getMoment("DD",  { id, id }, factorIds));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, 4. * popIdPower));

              if(popIdPower > 1)
              {
                factorIds = (*it)->getFactorIndices();
                sslib.dropFactorIds(factorIds, id, 1);

                col = sslib.findCompressedIndex(sslib.getMoment("Dr",  { id, id }, factorIds));
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, (popIdPower * (popIdPower - 1)) / 2.));
              }
            }
          }

          else // D_id_r_*
          {
            if(popIdPower > 1)
            {
              std::vector<size_t> factorIds = (*it)->getFactorIndices();
              sslib.dropFactorIds(factorIds, id, 2);

              col = sslib.findCompressedIndex(sslib.getMoment("Dr",  { id, id }, factorIds));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, (popIdPower * (popIdPower - 1)) / 2.));
            }
          }
        }
      }

      else if((*it)->getPrefix() == "DD")
      {
        int x = computeDDMainDiagContribution_(*it, id);
        coeffs.emplace_back(Eigen::Triplet<double>(row, row, x));

        if(popIdCount == 2)
        {
          std::vector<size_t> factorIds = (*it)->getFactorIndices();
          col = sslib.findCompressedIndex(sslib.getMoment("pi2", { id, id, id, id }, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));

          factorIds.push_back(id);
          col = sslib.findCompressedIndex(sslib.getMoment("Dr",  { id, id }, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));
        }

        if(popIdPower > 1)
        {
          std::vector<size_t> factorIds = (*it)->getFactorIndices();
          sslib.dropFactorIds(factorIds, id, 2);

          col = sslib.findCompressedIndex(sslib.getMoment("DD",  { id, id }, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, (popIdPower * (popIdPower - 1)) / 2.));
        }
      }

      else if((*it)->getPrefix() == "Hl")
      {
        if(popIdCount == 2)
        {
          coeffs.emplace_back(Eigen::Triplet<double>(row, row, (popIdPower * (popIdPower - 1)) / 2.));

          if(popIdPower > 1)
          {
            std::vector<size_t> factorIds = (*it)->getFactorIndices();
            sslib.dropFactorIds(factorIds, id, 2);

            col = sslib.findCompressedIndex(sslib.getMoment("Hl", { id, id }, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, (popIdPower * (popIdPower - 1)) / 2.));
          }
        }
      }

      else if((*it)->getPrefix() == "Hr")
      {
        if(popIdCount == 2)
          coeffs.emplace_back(Eigen::Triplet<double>(row, row, -1.));
      }

      else if((*it)->getPrefix() == "pi2")
      {
        int x = computePi2MainDiagContribution_(*it, id);
        coeffs.emplace_back(Eigen::Triplet<double>(row, row, x));

        std::vector<size_t> factorIds = (*it)->getFactorIndices();
        factorIds.push_back(id);

        col = sslib.findCompressedIndex(sslib.getMoment("Dr",  { id, id }, factorIds));
        coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1. + popIdPower/2.));

        if(popIdPower > 0)
        {
          sslib.dropFactorIds(factorIds, id, 2);
          col = sslib.findCompressedIndex(sslib.getMoment("Dr",  { id, id }, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, (-2. * popIdPower) / 4.));

          if(popIdPower > 1)
          {
            sslib.dropFactorIds(factorIds, id, 1);

            col = sslib.findCompressedIndex(sslib.getMoment("pi2", { id, id, id, id }, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, (popIdPower * (popIdPower - 1)) / 2.));
          }
        }
      }

      else if((*it)->getPrefix() != "I")
        throw bpp::Exception("Drift::mis-specified Moment prefix: " + (*it)->getPrefix());
    }

    Eigen::SparseMatrix<double> mat(sizeOfBasis, sizeOfBasis);
    mat.setFromTriplets(std::begin(coeffs), std::end(coeffs));
    mat.makeCompressed();
    mat *= (getParameterValue("1/2N_" + bpp::TextTools::toString(id)));
    matrices_.emplace_back(mat);
  }

  setIdentity_(sizeOfBasis);
  assembleTransitionMatrix_();
}

void Drift::updateMatrices_()
{
  for(size_t i = 0; i < matrices_.size(); ++i)
  {
    size_t id = popIndices_[i];
    std::string paramName = "1/2N_" + bpp::TextTools::toString(id);

    double prevVal = prevParams_.getParameterValue(paramName);
    double newVal = getParameterValue(paramName);

    if(newVal != prevVal)
      matrices_[i] *= (newVal / prevVal);
  }

  assembleTransitionMatrix_();
  prevParams_.matchParametersValues(getParameters());
}
