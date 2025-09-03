/*
 * Authors: Gustavo V. Barroso
 * Created: 09/08/2022
 * Last modified: 03/09/2025
 *
 */

#include <Bpp/Text/TextTools.h>

#include "Drift.hpp"

int Drift::computeDMainDiagContribution_(std::shared_ptr<Moment> mom, size_t id)
{
  size_t popIdCount = mom->countInstances(id);
  size_t popIdPower = mom->getPopFactorPower(id);
  size_t totalPopIdCount = popIdCount + popIdPower;

  int a = 0;
  int b = -1;

  for(size_t i = 0; i < totalPopIdCount + popIdCount - 1; ++i)
  {
    a += b;
    --b;
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
    if(mom->getPopIndices()[0] == id)
    {
      int a = 0;
      int b = -1;

      for(size_t i = 0; i < totalPopIdCount; ++i)
      {
        a += b;
        --b;
      }

      return a;
    }

    else
    {
      if(totalPopIdCount > 2)
      {
        int a = 0;
        int b = -1;

        for(size_t i = 0; i < totalPopIdCount - 2; ++i)
        {
          a += b;
          --b;
        }

        return a;
      }

      else
        return 0;
    }
  }

  else // if popIdCount == 0
  {
    if(totalPopIdCount > 1)
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

    else
      return 0;
  }
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
    int a = 0;
    int b = -1;

    for(size_t i = 0; i < totalPopIdCount; ++i)
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
      int a = 0;
      int b = -1;

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

  else if(popIdCount == 0 && popIdPower > 1)
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

  else
    return 0;
}

void Drift::setUpMatrices_(const SumStatsLibrary& sslib)
{
  size_t numPops = getParameters().size();
  size_t sizeOfBasis = sslib.getSizeOfBasis();
  matrices_.reserve(numPops);

  for(size_t i = 0; i < numPops; ++i)
  {
    size_t id = popIndices_[i];
    std::vector<Eigen::Triplet<mpfr::mpreal>> coeffs(0);
    coeffs.reserve(sizeOfBasis);

    for(auto it = std::begin(sslib.getBasis()); it != std::end(sslib.getBasis()); ++it)
    {
      int row = it - std::begin(sslib.getBasis());
      int col = -1;

      int popIdCount = static_cast<int>((*it)->countInstances(id));
      int popIdPower = static_cast<int>((*it)->getPopFactorPower(id));

      if((*it)->getPrefix() == "Dr")
      {
        int md = computeDrMainDiagContribution_(*it, id);
        coeffs.emplace_back(Eigen::Triplet<mpfr::mpreal>(row, row, md));

        if(popIdCount == 2)
        {
          if(popIdPower > 0)
          {
            std::vector<size_t> factorIds = (*it)->getFactorIndices();
            sslib.dropFactorIds(factorIds, id, 1);

            while(factorIds.size() > sslib.getFactorOrder()) // NOTE truncation
              factorIds.pop_back();

            col = sslib.findCompressedIndex(sslib.getMoment("DD",  { id, id }, factorIds));
            coeffs.emplace_back(Eigen::Triplet<mpfr::mpreal>(row, col, 4. * popIdPower));

            if(popIdPower > 1)
            {
              factorIds = (*it)->getFactorIndices();
              sslib.dropFactorIds(factorIds, id, 2);

              col = sslib.findCompressedIndex(sslib.getMoment("Dr",  { id, id }, factorIds));
              coeffs.emplace_back(Eigen::Triplet<mpfr::mpreal>(row, col, (popIdPower * (popIdPower - 1)) / 2.));
            }
          }
        }

        else if(popIdCount == 1) // D_id_r_*
        {
          if((*it)->getPopIndices()[0] == id && popIdPower > 1)
          {
            std::vector<size_t> factorIds = (*it)->getFactorIndices();
            sslib.dropFactorIds(factorIds, id, 2);

            col = sslib.findCompressedIndex(sslib.getMoment("Dr",  { id, sslib.fetchOtherId(id) }, factorIds));
            coeffs.emplace_back(Eigen::Triplet<mpfr::mpreal>(row, col, (popIdPower * (popIdPower - 1)) / 2.));
          }

          else if((*it)->getPopIndices()[0] != id && popIdPower > 0)
          {
            std::vector<size_t> factorIds = (*it)->getFactorIndices();
            sslib.dropFactorIds(factorIds, id, 1);

            while(factorIds.size() > sslib.getFactorOrder()) // NOTE truncation
              factorIds.pop_back();

            col = sslib.findCompressedIndex(sslib.getMoment("DD",  { id, sslib.fetchOtherId(id) }, factorIds));
            coeffs.emplace_back(Eigen::Triplet<mpfr::mpreal>(row, col, 4. * popIdPower));

            if(popIdPower > 1)
            {
              factorIds = (*it)->getFactorIndices();
              sslib.dropFactorIds(factorIds, id, 2);

              col = sslib.findCompressedIndex(sslib.getMoment("Dr",  { id, sslib.fetchOtherId(id) }, factorIds));
              coeffs.emplace_back(Eigen::Triplet<mpfr::mpreal>(row, col, (popIdPower * (popIdPower - 1)) / 2.));
            }
          }
        }

        else if(popIdCount == 0 && popIdPower > 1)
        {
          std::vector<size_t> factorIds = (*it)->getFactorIndices();
          sslib.dropFactorIds(factorIds, id, 2);

          col = sslib.findCompressedIndex(sslib.getMoment("Dr",  { sslib.fetchOtherId(id), sslib.fetchOtherId(id) }, factorIds));
          coeffs.emplace_back(Eigen::Triplet<mpfr::mpreal>(row, col, (popIdPower * (popIdPower - 1)) / 2.));
        }
      }

      #ifdef NAKED_D
      else if((*it)->getPrefix() == "D")
      {
        int md = computeDMainDiagContribution_(*it, id);
        coeffs.emplace_back(Eigen::Triplet<mpfr::mpreal>(row, row, md));

        if(popIdPower > 1)
        {
          std::vector<size_t> factorIds = (*it)->getFactorIndices();
          sslib.dropFactorIds(factorIds, id, 2);

          col = sslib.findCompressedIndex(sslib.getMoment("D", (*it)->getPopIndices(), factorIds));
          coeffs.emplace_back(Eigen::Triplet<mpfr::mpreal>(row, col, (popIdPower * (popIdPower - 1)) / 2.));
        }
      }
      #endif

      else if((*it)->getPrefix() == "DD")
      {
        int md = computeDDMainDiagContribution_(*it, id);
        coeffs.emplace_back(Eigen::Triplet<mpfr::mpreal>(row, row, md));

        if(popIdCount == 2)
        {
          std::vector<size_t> factorIds = (*it)->getFactorIndices();
          col = sslib.findCompressedIndex(sslib.getMoment("pi2", { id, id, id, id }, factorIds));
          coeffs.emplace_back(Eigen::Triplet<mpfr::mpreal>(row, col, 1.));

          factorIds.push_back(id);
          col = sslib.findCompressedIndex(sslib.getMoment("Dr",  { id, id }, factorIds));
          coeffs.emplace_back(Eigen::Triplet<mpfr::mpreal>(row, col, 1.));
        }

        if(popIdPower > 1)
        {
          std::vector<size_t> factorIds = (*it)->getFactorIndices();
          sslib.dropFactorIds(factorIds, id, 2);

          col = sslib.findCompressedIndex(sslib.getMoment("DD",  { id, id }, factorIds));
          coeffs.emplace_back(Eigen::Triplet<mpfr::mpreal>(row, col, (popIdPower * (popIdPower - 1)) / 2.));
        }
      }

      else if((*it)->getPrefix() == "Hl")
      {
        int totalPopIdCount = static_cast<int>(popIdCount) + static_cast<int>(popIdPower);
        coeffs.emplace_back(Eigen::Triplet<mpfr::mpreal>(row, row, -(totalPopIdCount * (totalPopIdCount - 1)) / 2.));

        if(popIdCount == 1 && popIdPower > 0)
        {
          int sign = std::pow(-1, (*it)->getPopIndices()[0] == id); // sign of contributions that would cancel out if p1(1-p0) == p0(1-p1)

          std::vector<size_t> factorIds = (*it)->getFactorIndices();
          sslib.dropFactorIds(factorIds, id, 1);

          col = sslib.findCompressedIndex(sslib.getMoment("Hl", (*it)->getPopIndices(), factorIds));
          coeffs.emplace_back(Eigen::Triplet<mpfr::mpreal>(row, col, sign * popIdPower));
        }

        if(popIdPower > 1)
        {
          std::vector<size_t> factorIds = (*it)->getFactorIndices();
          sslib.dropFactorIds(factorIds, id, 2);

          col = sslib.findCompressedIndex(sslib.getMoment("Hl", (*it)->getPopIndices(), factorIds));
          coeffs.emplace_back(Eigen::Triplet<mpfr::mpreal>(row, col, (popIdPower * (popIdPower - 1)) / 2.));
        }
      }

      else if((*it)->getPrefix() == "Hr")
      {
        if(popIdCount == 2)
          coeffs.emplace_back(Eigen::Triplet<mpfr::mpreal>(row, row, -1.));
      }

      else if((*it)->getPrefix() == "pi2")
      {
        auto tmpPi2 = std::dynamic_pointer_cast<Pi2Moment>(*it);
        assert(tmpPi2 != nullptr);

        size_t countLeft = tmpPi2->getLeftHetStat()->countInstances(id);
        size_t countRight = tmpPi2->getRightHetStat()->countInstances(id);

        // folding factor... WARNING only in D & Dr contributions to pi2 stats?! NOTE compare to moments.LD approach
        int f = static_cast<int>((*it)->getNumberOfAliases() + 1);

        int md = computePi2MainDiagContribution_(*it, id);
        coeffs.emplace_back(Eigen::Triplet<mpfr::mpreal>(row, row, md));

        if(countLeft == 2 && countRight == 2)
        {
          std::vector<size_t> factorIds = (*it)->getFactorIndices();
          factorIds.push_back(id);

          col = sslib.findCompressedIndex(sslib.getMoment("Dr",  { id, id }, factorIds));
          coeffs.emplace_back(Eigen::Triplet<mpfr::mpreal>(row, col, 1. + popIdPower / 2.));

          if(popIdPower > 0)
          {
            sslib.dropFactorIds(factorIds, id, 2);

            col = sslib.findCompressedIndex(sslib.getMoment("Dr",  { id, id }, factorIds));
            coeffs.emplace_back(Eigen::Triplet<mpfr::mpreal>(row, col, -popIdPower / 4.));

            if(popIdPower > 1)
            {
              col = sslib.findCompressedIndex(sslib.getMoment("pi2", (*it)->getPopIndices(), factorIds));
              coeffs.emplace_back(Eigen::Triplet<mpfr::mpreal>(row, col, (popIdPower * (popIdPower - 1)) / 2.));
            }
          }
        }

        else if(countLeft == 2 && countRight == 1)
        {
          std::vector<size_t> factorIds = (*it)->getFactorIndices();
          factorIds.push_back(id);

          col = sslib.findCompressedIndex(sslib.getMoment("Dr", { id, sslib.fetchOtherId(id) }, factorIds));
          coeffs.emplace_back(Eigen::Triplet<mpfr::mpreal>(row, col, f * (1. / 4. + popIdPower / 8.)));

          #ifdef NAKED_D
          int signD = std::pow(-1, (*it)->getPopIndices()[3] == id);
          col = sslib.findCompressedIndex(sslib.getMoment("D", { id }, factorIds));
          coeffs.emplace_back(Eigen::Triplet<mpfr::mpreal>(row, col, f * signD * (1. / 4. + popIdPower / 8.)));
          #endif

          if(popIdPower > 0)
          {
            sslib.dropFactorIds(factorIds, id, 2);

            col = sslib.findCompressedIndex(sslib.getMoment("Dr", { id, sslib.fetchOtherId(id) }, factorIds));
            coeffs.emplace_back(Eigen::Triplet<mpfr::mpreal>(row, col, -f * popIdPower / 8.));

            #ifdef NAKED_D
            col = sslib.findCompressedIndex(sslib.getMoment("D", { id }, factorIds));
            coeffs.emplace_back(Eigen::Triplet<mpfr::mpreal>(row, col, -f * signD * popIdPower / 8.));
            #endif

            if(popIdPower > 1)
            {
              sslib.dropFactorIds(factorIds, id, 1);

              col = sslib.findCompressedIndex(sslib.getMoment("pi2", (*it)->getPopIndices(), factorIds));
              coeffs.emplace_back(Eigen::Triplet<mpfr::mpreal>(row, col, (popIdPower * (popIdPower - 1)) / 2.));
            }
          }
        }

        else if(countLeft == 2 && countRight == 0)
        {
          if(popIdPower > 1)
          {
            std::vector<size_t> factorIds = (*it)->getFactorIndices();
            sslib.dropFactorIds(factorIds, id, 2);

            col = sslib.findCompressedIndex(sslib.getMoment("pi2", (*it)->getPopIndices(), factorIds));
            coeffs.emplace_back(Eigen::Triplet<mpfr::mpreal>(row, col, (popIdPower * (popIdPower - 1)) / 2.));
          }
        }

        else if(countLeft == 1 && countRight == 2) // NOTE check right-locus permutations when we include migration
        {
          int sign = std::pow(-1, (*it)->getPopIndices()[0] != id); // sign of contributions that would cancel out if p1(1-p0) == p0(1-p1)
          std::vector<size_t> factorIds = (*it)->getFactorIndices();

          col = sslib.findCompressedIndex(sslib.getMoment("Dr", { id, id }, factorIds));
          coeffs.emplace_back(Eigen::Triplet<mpfr::mpreal>(row, col, sign * (1. / 4. + popIdPower / 4.)));

          factorIds.push_back(sslib.fetchOtherId(id));

          col = sslib.findCompressedIndex(sslib.getMoment("Dr", { id, id }, factorIds));
          coeffs.emplace_back(Eigen::Triplet<mpfr::mpreal>(row, col, 1. / 4. + popIdPower / 4.));

          factorIds.pop_back();

          if(popIdPower > 0)
          {
            sslib.dropFactorIds(factorIds, id, 1);

            col = sslib.findCompressedIndex(sslib.getMoment("pi2", (*it)->getPopIndices(), factorIds));
            coeffs.emplace_back(Eigen::Triplet<mpfr::mpreal>(row, col, -sign * popIdPower));

            col = sslib.findCompressedIndex(sslib.getMoment("Dr", { id, id }, factorIds));
            coeffs.emplace_back(Eigen::Triplet<mpfr::mpreal>(row, col, -popIdPower / 4.));

            factorIds.push_back(sslib.fetchOtherId(id));

            col = sslib.findCompressedIndex(sslib.getMoment("Dr", { id, id }, factorIds));
            coeffs.emplace_back(Eigen::Triplet<mpfr::mpreal>(row, col, -sign * popIdPower / 4.));

            if(popIdPower > 1)
            {
              sslib.dropFactorIds(factorIds, id, 1);

              col = sslib.findCompressedIndex(sslib.getMoment("pi2", (*it)->getPopIndices(), factorIds));
              coeffs.emplace_back(Eigen::Triplet<mpfr::mpreal>(row, col, (popIdPower * (popIdPower - 1)) / 2.));
            }
          }
        }

        else if(countLeft == 1 && countRight == 1)
        {
          int sign = std::pow(-1, (*it)->getPopIndices()[0] != id); // sign of contributions that would cancel out if p1(1-p0) == p0(1-p1)
          std::vector<size_t> factorIds = (*it)->getFactorIndices();

          // (1-2p_id)^k
          col = sslib.findCompressedIndex(sslib.getMoment("Dr", { id, sslib.fetchOtherId(id) }, factorIds));
          coeffs.emplace_back(Eigen::Triplet<mpfr::mpreal>(row, col, f * sign * (popIdPower + 1) / 8.));

          /* cancel out due to right-locus permutation
          #ifdef NAKED_D
          int signD = std::pow(-1, (*it)->getPopIndices()[3] == id);
          col = sslib.findCompressedIndex(sslib.getMoment("D", { id }, factorIds));
          coeffs.emplace_back(Eigen::Triplet<mpfr::mpreal>(row, col, signD * (popIdPower + 1) / 8.));
          #endif
          */

          factorIds.push_back(sslib.fetchOtherId(id));

          col = sslib.findCompressedIndex(sslib.getMoment("Dr", { id, sslib.fetchOtherId(id) }, factorIds));
          coeffs.emplace_back(Eigen::Triplet<mpfr::mpreal>(row, col, f * (popIdPower + 1) / 8.));

          /* cancel out due to right-locus permutation
          #ifdef NAKED_D
          col = sslib.findCompressedIndex(sslib.getMoment("D", { id }, factorIds));
          coeffs.emplace_back(Eigen::Triplet<mpfr::mpreal>(row, col, signD * (popIdPower + 1) / 8.));
          #endif
          */

          factorIds.pop_back();

          if(popIdPower > 0)
          {
            // (1-2p_id)^(k-1)
            sslib.dropFactorIds(factorIds, id, 1);

            col = sslib.findCompressedIndex(sslib.getMoment("Dr", { id, sslib.fetchOtherId(id) }, factorIds));
            coeffs.emplace_back(Eigen::Triplet<mpfr::mpreal>(row, col, f * sign * (popIdPower + 1) / 8.));

            /* cancel out due to right-locus permutation
            #ifdef NAKED_D
            col = sslib.findCompressedIndex(sslib.getMoment("D", { id }, factorIds));
            coeffs.emplace_back(Eigen::Triplet<mpfr::mpreal>(row, col, sign * (popIdPower + 1) / 8.));
            #endif
            */

            factorIds.push_back(sslib.fetchOtherId(id));

            col = sslib.findCompressedIndex(sslib.getMoment("Dr", { id, sslib.fetchOtherId(id) }, factorIds));
            coeffs.emplace_back(Eigen::Triplet<mpfr::mpreal>(row, col, f * sign * (popIdPower - 1) / 8.));

            /* cancel out due to right-locus permutation
            #ifdef NAKED_D
            col = sslib.findCompressedIndex(sslib.getMoment("D", { id }, factorIds));
            coeffs.emplace_back(Eigen::Triplet<mpfr::mpreal>(row, col,  signD * (popIdPower + 1) / 8.));
            #endif
            */

            factorIds.pop_back();

            col = sslib.findCompressedIndex(sslib.getMoment("pi2", (*it)->getPopIndices(), factorIds));
            coeffs.emplace_back(Eigen::Triplet<mpfr::mpreal>(row, col, -sign * popIdPower));

            if(popIdPower > 1)
            {
              sslib.dropFactorIds(factorIds, id, 1);

              col = sslib.findCompressedIndex(sslib.getMoment("pi2", (*it)->getPopIndices(), factorIds));
              coeffs.emplace_back(Eigen::Triplet<mpfr::mpreal>(row, col, (popIdPower * (popIdPower - 1)) / 2.));
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
            col = sslib.findCompressedIndex(sslib.getMoment("pi2", (*it)->getPopIndices(), factorIds));
            coeffs.emplace_back(Eigen::Triplet<mpfr::mpreal>(row, col, sign * popIdPower));

            if(popIdPower > 1)
            {
              sslib.dropFactorIds(factorIds, id, 1);

              col = sslib.findCompressedIndex(sslib.getMoment("pi2", (*it)->getPopIndices(), factorIds));
              coeffs.emplace_back(Eigen::Triplet<mpfr::mpreal>(row, col, (popIdPower * (popIdPower - 1)) / 2.));
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
            coeffs.emplace_back(Eigen::Triplet<mpfr::mpreal>(row, col, -popIdPower / 4.));

            factorIds.push_back(sslib.fetchOtherId(id));
            factorIds.push_back(sslib.fetchOtherId(id));

            col = sslib.findCompressedIndex(sslib.getMoment("Dr", { id, id }, factorIds));
            coeffs.emplace_back(Eigen::Triplet<mpfr::mpreal>(row, col, popIdPower / 4.));

            if(popIdPower > 1)
            {
              sslib.dropFactorIds(factorIds, sslib.fetchOtherId(id), 2);
              sslib.dropFactorIds(factorIds, id, 1);

              col = sslib.findCompressedIndex(sslib.getMoment("pi2", (*it)->getPopIndices(), factorIds));
              coeffs.emplace_back(Eigen::Triplet<mpfr::mpreal>(row, col, (popIdPower * (popIdPower - 1)) / 2.));
            }
          }
        }

        else if(countLeft == 0 && countRight == 1)
        {
          if(popIdPower > 0)
          {
            std::vector<size_t> factorIds = (*it)->getFactorIndices();
            sslib.dropFactorIds(factorIds, id, 1);

            col = sslib.findCompressedIndex(sslib.getMoment("Dr", { id, sslib.fetchOtherId(id) }, factorIds));
            coeffs.emplace_back(Eigen::Triplet<mpfr::mpreal>(row, col, -f * popIdPower / 8.));

            #ifdef NAKED_D
            int signD = std::pow(-1, (*it)->getPopIndices()[3] == id);
            col = sslib.findCompressedIndex(sslib.getMoment("D", { id }, factorIds));
            coeffs.emplace_back(Eigen::Triplet<mpfr::mpreal>(row, col, -f * signD * (-popIdPower) / 8.));
            #endif

            factorIds.push_back(sslib.fetchOtherId(id));
            factorIds.push_back(sslib.fetchOtherId(id));

            col = sslib.findCompressedIndex(sslib.getMoment("Dr", { id, sslib.fetchOtherId(id) }, factorIds));
            coeffs.emplace_back(Eigen::Triplet<mpfr::mpreal>(row, col, f * popIdPower / 8.));

            #ifdef NAKED_D
            col = sslib.findCompressedIndex(sslib.getMoment("D", { id }, factorIds));
            coeffs.emplace_back(Eigen::Triplet<mpfr::mpreal>(row, col, f * signD * (-popIdPower) / 8.));
            #endif

            if(popIdPower > 1)
            {
              sslib.dropFactorIds(factorIds, id, 1);
              sslib.dropFactorIds(factorIds, sslib.fetchOtherId(id), 2);

              col = sslib.findCompressedIndex(sslib.getMoment("pi2", (*it)->getPopIndices(), factorIds));
              coeffs.emplace_back(Eigen::Triplet<mpfr::mpreal>(row, col, (popIdPower * (popIdPower - 1)) / 2.));
            }
          }
        }

        else if(popIdCount == 0 && popIdPower > 1)
        {
          std::vector<size_t> factorIds = (*it)->getFactorIndices();
          sslib.dropFactorIds(factorIds, id, 2);

          col = sslib.findCompressedIndex(sslib.getMoment("pi2", (*it)->getPopIndices(), factorIds));
          coeffs.emplace_back(Eigen::Triplet<mpfr::mpreal>(row, col, (popIdPower * (popIdPower - 1)) / 2.));
        }
      }

      else if((*it)->getPrefix() != "I")
        throw bpp::Exception("Drift::mis-specified Moment prefix: " + (*it)->getPrefix());
    }

    Eigen::SparseMatrix<mpfr::mpreal> mat(sizeOfBasis, sizeOfBasis);
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

    mpfr::mpreal prevVal = prevParams_.getParameterValue(paramName);
    mpfr::mpreal newVal = getParameterValue(paramName);

    if(newVal != prevVal)
      matrices_[i] *= (newVal / prevVal);
  }

  assembleTransitionMatrix_();
  prevParams_.matchParametersValues(getParameters());
}
