/*
 * Authors: Gustavo V. Barroso
 * Created: 09/08/2022
 * Last modified: 03/10/2023
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
  auto tmpPi2 = std::dynamic_pointer_cast<Pi2Moment>(mom);

  size_t countLeft = tmpPi2->getLeftHetStat()->countInstances(id);
  size_t countRight = tmpPi2->getRightHetStat()->countInstances(id);
  size_t popIdCount = countLeft + countRight;
  size_t popIdPower = tmpPi2->getPopFactorPower(id);
  size_t totalPopIdCount = popIdCount + popIdPower;

  if(popIdCount == 4)
  {
    int a = -1;
    int b = -1;

    for(size_t i = 0; i < totalPopIdCount - 3; ++i)
    {
      a += b;
      --b;
    }

    return a;
  }

  else if((countLeft == 2) || (countRight == 2))
  {
    if(popIdCount == 3)
    {
      int a = -1;
      int b = 0;

      for(size_t i = 0; i < totalPopIdCount - 2; ++i)
      {
        a += b;
        --b;
      }

      return a;
    }

    else
    {
      int a = -1;
      int b = 0;

      for(size_t i = 0; i < totalPopIdCount - popIdCount; ++i)
      {
        a += b;
        --b;
      }

      return a;
    }
  }

  else if((countLeft == 1) || (countRight == 1))
  {
    int a = 0;
    int b = -1;

    for(size_t i = 0; i < totalPopIdCount - popIdCount; ++i)
    {
      a += b;
      --b;
    }

    return a;
  }

  else
    return 0;
}

int Drift::computePi2OffDiagContribution_(std::shared_ptr<Moment> mom, size_t id)
{
  size_t totalPopIdCount = mom->countInstances(id) + mom->getPopFactorPower(id);

  int a = 0;
  int b = 1;

  for(size_t i = 0; i < totalPopIdCount - 2; ++i)
  {
    a += b;
    ++b;
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

      int popIdCount = static_cast<int>((*it)->countInstances(id));
      int popIdPower = static_cast<int>((*it)->getPopFactorPower(id));

      (*it)->printAttributes(std::cout);

      if((*it)->getPrefix() == "D")
      {
        if(popIdCount == 1)
        {
          int md = computeDMainDiagContribution_(*it, id);
          coeffs.emplace_back(Eigen::Triplet<double>(row, row, md));

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
          int md = computeDrMainDiagContribution_(*it, id);
          coeffs.emplace_back(Eigen::Triplet<double>(row, row, md));

          if(popIdCount == 2) // D_id_r_id
          {
            if(popIdPower > 0)
            {
              std::vector<size_t> factorIds = (*it)->getFactorIndices();
              sslib.dropFactorIds(factorIds, id, 1);

              if((popIdPower > static_cast<int>((sslib.getFactorOrder() + 1))) || (((*it)->getFactorPower() > sslib.getFactorOrder()) && popIdPower > 1))
                sslib.dropFactorIds(factorIds, id, 1);

              while(factorIds.size() > sslib.getFactorOrder()) // NOTE
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
        int md = computeDDMainDiagContribution_(*it, id);
        coeffs.emplace_back(Eigen::Triplet<double>(row, row, md));

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
        auto tmpPi2 = std::dynamic_pointer_cast<Pi2Moment>(*it);
        assert(tmpPi2 != nullptr);

        size_t countLeft = tmpPi2->getLeftHetStat()->countInstances(id);
        size_t countRight = tmpPi2->getRightHetStat()->countInstances(id);

        int md = computePi2MainDiagContribution_(*it, id);
        coeffs.emplace_back(Eigen::Triplet<double>(row, row, md));

        if(countLeft == 2 && countRight == 2)
        {
          std::vector<size_t> factorIds = (*it)->getFactorIndices();
          factorIds.push_back(id);

          col = sslib.findCompressedIndex(sslib.getMoment("Dr",  { id, id }, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1. + popIdPower / 2.));

          if(popIdPower > 0)
          {
            sslib.dropFactorIds(factorIds, id, 1);
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

        else if(countLeft == 2 && countRight == 1)
        {
          std::vector<size_t> factorIds = (*it)->getFactorIndices();
          factorIds.push_back(id);

          col = sslib.findCompressedIndex(sslib.getMoment("Dr", { id, sslib.fetchOtherId(id) }, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1. / 4. + popIdPower / 8.));

          // NOTE left out because of right-locus symmetries
          //col = sslib.findCompressedIndex(sslib.getMoment("D", { id }, factorIds));
          //coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1. / 4. + popIdPower / 8.));

          if(popIdPower > 0)
          {
            sslib.dropFactorIds(factorIds, id, 2);

            col = sslib.findCompressedIndex(sslib.getMoment("Dr", { id, sslib.fetchOtherId(id) }, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, -popIdPower / 8.));

            // NOTE left out because of right-locus symmetries
            //col = sslib.findCompressedIndex(sslib.getMoment("D", { id }, factorIds));
            //coeffs.emplace_back(Eigen::Triplet<double>(row, col, -popIdPower / 8.));

            if(popIdPower > 1)
            {
              sslib.dropFactorIds(factorIds, id, 1);

              int y = computePi2MainDiagContribution_(sslib.getMoment("pi2", (*it)->getPopIndices(), factorIds), id); // NOTE
              col = sslib.findCompressedIndex(sslib.getMoment("pi2", (*it)->getPopIndices(), factorIds));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, y));
            }
          }
        }

        else if(countLeft == 2 && countRight == 0)
        {
          if(popIdPower > 1)
          {
            std::vector<size_t> factorIds = (*it)->getFactorIndices();
            sslib.dropFactorIds(factorIds, id, 2);

            int x = computePi2OffDiagContribution_(sslib.getMoment("pi2", (*it)->getPopIndices(), factorIds), id);
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, x));
          }
        }

        else if(countLeft == 1 && countRight == 2)
        {
          int sign = std::pow(-1, (*it)->getPopIndices()[0] != id); // sign of contributions that would cancel out if p1(1-p0) == p0(1-p1)
          std::vector<size_t> factorIds = (*it)->getFactorIndices();

          col = sslib.findCompressedIndex(sslib.getMoment("Dr", { id, id }, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, sign * (1. / 4. + popIdPower / 4.)));

          factorIds.push_back(sslib.fetchOtherId(id));

          col = sslib.findCompressedIndex(sslib.getMoment("Dr", { id, id }, factorIds));
          coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1. / 4. + popIdPower / 4.));

          factorIds.pop_back();

          if(popIdPower > 0)
          {
            sslib.dropFactorIds(factorIds, id, 1);

            col = sslib.findCompressedIndex(sslib.getMoment("pi2", (*it)->getPopIndices(), factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, -sign * popIdPower));

            col = sslib.findCompressedIndex(sslib.getMoment("Dr", { id, id }, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1. / 4. + popIdPower / 4.));

            factorIds.push_back(sslib.fetchOtherId(id));

            col = sslib.findCompressedIndex(sslib.getMoment("Dr", { id, id }, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, sign * (1. / 4. + popIdPower / 4.)));

            if(popIdPower > 1)
            {
              sslib.dropFactorIds(factorIds, id, 1);

              int x = computePi2OffDiagContribution_(sslib.getMoment("pi2", (*it)->getPopIndices(), factorIds), id);
              col = sslib.findCompressedIndex(sslib.getMoment("pi2", (*it)->getPopIndices(), factorIds));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, x));
            }
          }
        }

        else if(countLeft == 1 && countRight == 1) // a bit convoluted because popIdPower == 0 is a weird special case!
        {
          if(popIdPower == 0)
          {
            int sign = std::pow(-1, popIdPower); // sign of contributions that would cancel out if p1(1-p0) == p0(1-p1)
            std::vector<size_t> factorIds = (*it)->getFactorIndices();

            col = sslib.findCompressedIndex(sslib.getMoment("Dr", { id, sslib.fetchOtherId(id) }, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, sign * (1. / 8. + popIdPower / 8.)));

            col = sslib.findCompressedIndex(sslib.getMoment("D", { id }, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, sign * (1. / 8. + popIdPower / 8.)));

            factorIds.push_back(sslib.fetchOtherId(id));

            col = sslib.findCompressedIndex(sslib.getMoment("Dr", { id, sslib.fetchOtherId(id) }, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1. / 8. + popIdPower / 8.));

            col = sslib.findCompressedIndex(sslib.getMoment("D", { id }, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col,  sign * (1. / 8. + popIdPower / 8.)));
          }

          else if(popIdPower > 0)
          {
            std::vector<size_t> factorIds = (*it)->getFactorIndices();

            // (1-2p)^k
            col = sslib.findCompressedIndex(sslib.getMoment("Dr", { id, sslib.fetchOtherId(id) }, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, (popIdPower + 1) / 8.));

            col = sslib.findCompressedIndex(sslib.getMoment("D", { id }, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, (popIdPower + 1) / 8.));

            factorIds.push_back(sslib.fetchOtherId(id));

            col = sslib.findCompressedIndex(sslib.getMoment("Dr", { id, sslib.fetchOtherId(id) }, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, (popIdPower + 1) / 8.));

            col = sslib.findCompressedIndex(sslib.getMoment("D", { id }, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col,  (popIdPower + 1) / 8.));

            factorIds.pop_back();

            // (1-2p)^(k-1)
            sslib.dropFactorIds(factorIds, id, 1);

            col = sslib.findCompressedIndex(sslib.getMoment("Dr", { id, sslib.fetchOtherId(id) }, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, (popIdPower + 1) / 8.));

            col = sslib.findCompressedIndex(sslib.getMoment("D", { id }, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, (popIdPower + 1) / 8.));

            factorIds.push_back(sslib.fetchOtherId(id));

            col = sslib.findCompressedIndex(sslib.getMoment("Dr", { id, sslib.fetchOtherId(id) }, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, (popIdPower + 1) / 8.));

            col = sslib.findCompressedIndex(sslib.getMoment("D", { id }, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col,  (popIdPower + 1) / 8.));

            factorIds.pop_back();

            // other pi2
            col = sslib.findCompressedIndex(sslib.getMoment("pi2", (*it)->getPopIndices(), factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, -popIdPower));

            if(popIdPower > 1)
            {
              sslib.dropFactorIds(factorIds, id, 1);

              col = sslib.findCompressedIndex(sslib.getMoment("pi2", (*it)->getPopIndices(), factorIds));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, popIdPower));
            }
          }
        }

        else if(countLeft == 1 && countRight == 0)
        {
          if(popIdPower > 0)
          {
            std::vector<size_t> factorIds = (*it)->getFactorIndices();
            sslib.dropFactorIds(factorIds, id, 1);

            int sign = std::pow(-1, (*it)->getPopIndices()[0] == id); // sign of contributions that would cancel out if p1(1-p0) == p0(1-p1)
            col = sslib.findCompressedIndex(sslib.getMoment("pi2", (*it)->getPopIndices(), factorIds), id);
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, sign * popIdPower));

            if(popIdPower > 1)
            {
              std::vector<size_t> factorIds = (*it)->getFactorIndices();
              sslib.dropFactorIds(factorIds, id, 1);

              int x = computePi2OffDiagContribution_(sslib.getMoment("pi2", (*it)->getPopIndices(), factorIds), id);
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, x));
            }
          }
        }

        else if(countLeft == 0 && countRight == 2)
        {
          if(popIdPower > 0)
          {
            std::vector<size_t> factorIds = (*it)->getFactorIndices();
            sslib.dropFactorIds(factorIds, id, 1);

            col = sslib.findCompressedIndex(sslib.getMoment("Dr", { id, id }, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, -(1. / 4. + (popIdPower - 1) / 4.)));

            factorIds.push_back(sslib.fetchOtherId(id));
            factorIds.push_back(sslib.fetchOtherId(id));

            col = sslib.findCompressedIndex(sslib.getMoment("Dr", { id, id }, factorIds));
            coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1. / 4. + (popIdPower - 1) / 4.));

            if(popIdPower > 1)
            {
              sslib.dropFactorIds(factorIds, sslib.fetchOtherId(id), 2);
              sslib.dropFactorIds(factorIds, id, 1);

              int x = computePi2OffDiagContribution_(sslib.getMoment("pi2", (*it)->getPopIndices(), factorIds), id);
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, x));
            }
          }
        }

        else if(countLeft == 0 && countRight == 1)
        {
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
