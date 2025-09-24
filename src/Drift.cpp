/*
 * Authors: Gustavo V. Barroso
 * Created: 09/08/2022
 * Last modified: 24/09/2025
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
  size_t basisSize = sslib.getSizeOfBasis();
  matrices_.reserve(numPops);

  const auto& basis = sslib.getBasis();
  const size_t numThreads = omp_get_max_threads();

  std::visit(overloaded{[&](auto const& mat)
  {
    using Scalar = typename std::decay_t<decltype(mat)>::Scalar;

    for(size_t i = 0; i < numPops; ++i)
    {
      const size_t id = popIndices_[i];
      const std::string paramName = "1/2N_" + bpp::TextTools::toString(id);
      const double coalRate = getParameterValue(paramName);

      std::vector<std::vector<Eigen::Triplet<Scalar>>> threadTriplets(numThreads);

      #pragma omp parallel for
      for(size_t row = 0; row < basisSize; ++row)
      {
        int col = -1;

        const auto& moment = basis[row];
        const std::string& prefix = moment->getPrefix();
        const size_t tid = omp_get_thread_num();
        auto& localTriplets = threadTriplets[tid];

        int popIdCount = static_cast<int>(moment->countInstances(id));
        int popIdPower = static_cast<int>(moment->getPopFactorPower(id));

        if(prefix == "Dr")
        {
          int md = computeDrMainDiagContribution_(moment, id);
          localTriplets.emplace_back(row, row, Scalar(md));

          if(popIdCount == 2)
          {
            if(popIdPower > 0)
            {
              std::vector<size_t> factorIds = moment->getFactorIndices();
              sslib.dropFactorIds(factorIds, id, 1);

              while(factorIds.size() > sslib.getFactorOrder()) // NOTE truncation
                factorIds.pop_back();

              col = sslib.findCompressedIndex(sslib.getMoment("DD", {id, id}, factorIds));
              localTriplets.emplace_back(row, col, Scalar(4. * popIdPower));

              if(popIdPower > 1)
              {
                factorIds = moment->getFactorIndices();
                sslib.dropFactorIds(factorIds, id, 2);

                col = sslib.findCompressedIndex(sslib.getMoment("Dr", {id, id}, factorIds));
                localTriplets.emplace_back(row, col, Scalar((popIdPower * (popIdPower - 1)) / 2.));
              }
            }
          }

          else if(popIdCount == 1) // D_id_r_*
          {
            if(moment->getPopIndices()[0] == id && popIdPower > 1)
            {
              std::vector<size_t> factorIds = moment->getFactorIndices();
              sslib.dropFactorIds(factorIds, id, 2);

              col = sslib.findCompressedIndex(sslib.getMoment("Dr", {id, sslib.fetchOtherId(id)}, factorIds));
              localTriplets.emplace_back(row, col, Scalar((popIdPower * (popIdPower - 1)) / 2.));
            }

            else if(moment->getPopIndices()[0] != id && popIdPower > 0)
            {
              std::vector<size_t> factorIds = moment->getFactorIndices();
              sslib.dropFactorIds(factorIds, id, 1);

              while(factorIds.size() > sslib.getFactorOrder()) // NOTE truncation
                factorIds.pop_back();

              col = sslib.findCompressedIndex(
              sslib.getMoment("DD", {id, sslib.fetchOtherId(id)}, factorIds));
              localTriplets.emplace_back(row, col, Scalar(4. * popIdPower));

              if(popIdPower > 1)
              {
                factorIds = moment->getFactorIndices();
                sslib.dropFactorIds(factorIds, id, 2);

                col = sslib.findCompressedIndex(
                sslib.getMoment("Dr", {id, sslib.fetchOtherId(id)}, factorIds));
                localTriplets.emplace_back(row, col, Scalar((popIdPower * (popIdPower - 1)) / 2.));
              }
            }
          }

          else if(popIdCount == 0 && popIdPower > 1)
          {
            std::vector<size_t> factorIds = moment->getFactorIndices();
            sslib.dropFactorIds(factorIds, id, 2);

            col = sslib.findCompressedIndex(sslib.getMoment("Dr", {sslib.fetchOtherId(id), sslib.fetchOtherId(id)}, factorIds));
            localTriplets.emplace_back(row, col, Scalar((popIdPower * (popIdPower - 1)) / 2.));
          }
        } // ends Dr prefix

        #ifdef NAKED_D
        else if(prefix == "D")
        {
          int md = computeDMainDiagContribution_(moment, id);
          localTriplets.emplace_back(row, row, Scalar(md));

          if(popIdPower > 1)
          {
            std::vector<size_t> factorIds = moment->getFactorIndices();
            sslib.dropFactorIds(factorIds, id, 2);

            col = sslib.findCompressedIndex(sslib.getMoment("D", moment->getPopIndices(), factorIds));
            localTriplets.emplace_back(row, col, Scalar((popIdPower * (popIdPower - 1)) / 2.));
          }
        }
        #endif

        else if(prefix == "DD")
        {
          int md = computeDDMainDiagContribution_(moment, id);
          localTriplets.emplace_back(row, row, Scalar(md));

          if(popIdCount == 2)
          {
            std::vector<size_t> factorIds = moment->getFactorIndices();
            col = sslib.findCompressedIndex(sslib.getMoment("pi2", {id, id, id, id}, factorIds));
            localTriplets.emplace_back(row, col, Scalar(1.));

            factorIds.push_back(id);
            col = sslib.findCompressedIndex(sslib.getMoment("Dr", {id, id}, factorIds));
            localTriplets.emplace_back(row, col, Scalar(1.));
          }

          if(popIdPower > 1)
          {
            std::vector<size_t> factorIds = moment->getFactorIndices();
            sslib.dropFactorIds(factorIds, id, 2);

            col = sslib.findCompressedIndex(sslib.getMoment("DD", {id, id}, factorIds));
            localTriplets.emplace_back(row, col, Scalar((popIdPower * (popIdPower - 1)) / 2.));
          }
        } // ends DD prefix

        else if(prefix == "Hl")
        {
          int totalPopIdCount = static_cast<int>(popIdCount) + static_cast<int>(popIdPower);
          localTriplets.emplace_back(row, row, Scalar(-(totalPopIdCount * (totalPopIdCount - 1)) / 2.));

          if(popIdCount == 1 && popIdPower > 0)
          {
            int sign = std::pow(-1, moment->getPopIndices()[0] == id); // sign of contributions that would cancel out if p1(1-p0) == p0(1-p1)

            std::vector<size_t> factorIds = moment->getFactorIndices();
            sslib.dropFactorIds(factorIds, id, 1);

            col = sslib.findCompressedIndex(sslib.getMoment("Hl", moment->getPopIndices(), factorIds));
            localTriplets.emplace_back(row, col, Scalar(sign * popIdPower));
          }

          if(popIdPower > 1)
          {
            std::vector<size_t> factorIds = moment->getFactorIndices();
            sslib.dropFactorIds(factorIds, id, 2);

            col = sslib.findCompressedIndex(sslib.getMoment("Hl", moment->getPopIndices(), factorIds));
            localTriplets.emplace_back(row, col, Scalar((popIdPower * (popIdPower - 1)) / 2.));
          }
        } // ends Hl prefix

        else if(prefix == "Hr")
        {
          if(popIdCount == 2)
            localTriplets.emplace_back(row, row, Scalar(-1.));
        } // ends Hr prefix

        else if(prefix == "pi2")
        {
          auto tmpPi2 = std::dynamic_pointer_cast<Pi2Moment>(moment);
          assert(tmpPi2 != nullptr);

          size_t countLeft = tmpPi2->getLeftHetStat()->countInstances(id);
          size_t countRight = tmpPi2->getRightHetStat()->countInstances(id);

          // folding factor... NOTE only in D & Dr contributions to pi2 stats?! --> compare to moments.LD approach
          int f = static_cast<int>(moment->getNumberOfAliases() + 1);

          int md = computePi2MainDiagContribution_(moment, id);
          localTriplets.emplace_back(row, row, Scalar(md));

          if(countLeft == 2 && countRight == 2)
          {
            std::vector<size_t> factorIds = moment->getFactorIndices();
            factorIds.push_back(id);

            col = sslib.findCompressedIndex(sslib.getMoment("Dr", {id, id}, factorIds));
            localTriplets.emplace_back(row, col, Scalar(1. + popIdPower / 2.));

            if(popIdPower > 0)
            {
              sslib.dropFactorIds(factorIds, id, 2);

              col = sslib.findCompressedIndex(sslib.getMoment("Dr", {id, id}, factorIds));
              localTriplets.emplace_back(row, col, Scalar(-popIdPower / 4.));

              if(popIdPower > 1)
              {
                col = sslib.findCompressedIndex(sslib.getMoment("pi2", moment->getPopIndices(), factorIds));
                localTriplets.emplace_back(row, col, Scalar((popIdPower * (popIdPower - 1)) / 2.));
              }
            }
          }

          else if(countLeft == 2 && countRight == 1)
          {
            std::vector<size_t> factorIds = moment->getFactorIndices();
            factorIds.push_back(id);

            col = sslib.findCompressedIndex(sslib.getMoment("Dr", {id, sslib.fetchOtherId(id)}, factorIds));
            localTriplets.emplace_back(row, col, f * Scalar((1. / 4. + popIdPower / 8.)));

            #ifdef NAKED_D
            int signD = std::pow(-1, moment->getPopIndices()[3] == id);
            col = sslib.findCompressedIndex(sslib.getMoment("D", {id}, factorIds));
            localTriplets.emplace_back(row, col, f * signD * Scalar((1. / 4. + popIdPower / 8.)));
            #endif

            if(popIdPower > 0)
            {
              sslib.dropFactorIds(factorIds, id, 2);

              col = sslib.findCompressedIndex(sslib.getMoment("Dr", {id, sslib.fetchOtherId(id)}, factorIds));
              localTriplets.emplace_back(row, col, Scalar(-f * popIdPower / 8.));

              #ifdef NAKED_D
              col = sslib.findCompressedIndex(sslib.getMoment("D", {id}, factorIds));
              localTriplets.emplace_back(row, col, Scalar(-f * signD * popIdPower / 8.));
              #endif

              if(popIdPower > 1)
              {
                sslib.dropFactorIds(factorIds, id, 1);

                col = sslib.findCompressedIndex(sslib.getMoment("pi2", moment->getPopIndices(), factorIds));
                localTriplets.emplace_back(row, col, Scalar((popIdPower * (popIdPower - 1)) / 2.));
              }
            }
          }

          else if(countLeft == 2 && countRight == 0)
          {
            if(popIdPower > 1)
            {
              std::vector<size_t> factorIds = moment->getFactorIndices();
              sslib.dropFactorIds(factorIds, id, 2);

              col = sslib.findCompressedIndex(sslib.getMoment("pi2", moment->getPopIndices(), factorIds));
              localTriplets.emplace_back(row, col, Scalar((popIdPower * (popIdPower - 1)) / 2.));
            }
          }

          else if(countLeft == 1 && countRight == 2) // NOTE check right-locus permutations when we include migration
          {
            int sign = std::pow(-1, moment->getPopIndices()[0] != id); // sign of contributions that would cancel out if p1(1-p0) == p0(1-p1)
            std::vector<size_t> factorIds = moment->getFactorIndices();

            col = sslib.findCompressedIndex(sslib.getMoment("Dr", {id, id}, factorIds));
            localTriplets.emplace_back(row, col, sign * Scalar((1. / 4. + popIdPower / 4.)));

            factorIds.push_back(sslib.fetchOtherId(id));

            col = sslib.findCompressedIndex(sslib.getMoment("Dr", {id, id}, factorIds));
            localTriplets.emplace_back(row, col, Scalar(1. / 4. + popIdPower / 4.));

            factorIds.pop_back();

            if(popIdPower > 0)
            {
              sslib.dropFactorIds(factorIds, id, 1);

              col = sslib.findCompressedIndex(sslib.getMoment("pi2", moment->getPopIndices(), factorIds));
              localTriplets.emplace_back(row, col, Scalar(-sign * popIdPower));

              col = sslib.findCompressedIndex(sslib.getMoment("Dr", {id, id}, factorIds));
              localTriplets.emplace_back(row, col, Scalar(-popIdPower / 4.));

              factorIds.push_back(sslib.fetchOtherId(id));

              col = sslib.findCompressedIndex(sslib.getMoment("Dr", {id, id}, factorIds));
              localTriplets.emplace_back(row, col, Scalar(-sign * popIdPower / 4.));

              if(popIdPower > 1)
              {
                sslib.dropFactorIds(factorIds, id, 1);

                col = sslib.findCompressedIndex(sslib.getMoment("pi2", moment->getPopIndices(), factorIds));
                localTriplets.emplace_back(row, col, Scalar((popIdPower * (popIdPower - 1)) / 2.));
              }
            }
          }

          else if(countLeft == 1 && countRight == 1)
          {
            int sign = std::pow(-1, moment->getPopIndices()[0] != id); // sign of contributions that would cancel out if p1(1-p0) == p0(1-p1)
            std::vector<size_t> factorIds = moment->getFactorIndices();

            // (1-2p_id)^k
            col = sslib.findCompressedIndex(sslib.getMoment("Dr", {id, sslib.fetchOtherId(id)}, factorIds));
            localTriplets.emplace_back(row, col, Scalar(f * sign * (popIdPower + 1) / 8.));

            /* cancel out due to right-locus permutation
            #ifdef NAKED_D
            int signD = std::pow(-1, moment->getPopIndices()[3] == id);
            col = sslib.findCompressedIndex(sslib.getMoment("D", { id }, factorIds));
            localTriplets.emplace_back(row, col, Scalar(signD * (popIdPower + 1) / 8.));
            #endif
            */

            factorIds.push_back(sslib.fetchOtherId(id));

            col = sslib.findCompressedIndex(
                sslib.getMoment("Dr", {id, sslib.fetchOtherId(id)}, factorIds));
            localTriplets.emplace_back(row, col, Scalar(f * (popIdPower + 1) / 8.));

            /* cancel out due to right-locus permutation
            #ifdef NAKED_D
            col = sslib.findCompressedIndex(sslib.getMoment("D", { id }, factorIds));
            localTriplets.emplace_back(row, col, Scalar(signD * (popIdPower + 1) / 8.));
            #endif
            */

            factorIds.pop_back();

            if(popIdPower > 0)
            {
              // (1-2p_id)^(k-1)
              sslib.dropFactorIds(factorIds, id, 1);

              col = sslib.findCompressedIndex(sslib.getMoment("Dr", {id, sslib.fetchOtherId(id)}, factorIds));
              localTriplets.emplace_back(row, col, Scalar(f * sign * (popIdPower + 1) / 8.));

              /* cancel out due to right-locus permutation
              #ifdef NAKED_D
              col = sslib.findCompressedIndex(sslib.getMoment("D", { id }, factorIds));
              localTriplets.emplace_back(row, col, Scalar(sign * (popIdPower + 1) / 8.));
              #endif
              */

              factorIds.push_back(sslib.fetchOtherId(id));

              col = sslib.findCompressedIndex(sslib.getMoment("Dr", {id, sslib.fetchOtherId(id)}, factorIds));
              localTriplets.emplace_back(row, col, Scalar(f * sign * (popIdPower - 1) / 8.));

              /* cancel out due to right-locus permutation
              #ifdef NAKED_D
              col = sslib.findCompressedIndex(sslib.getMoment("D", { id }, factorIds));
              localTriplets.emplace_back(row, col, Scalar(signD * (popIdPower + 1) / 8.));
              #endif
              */

              factorIds.pop_back();

              col = sslib.findCompressedIndex(sslib.getMoment("pi2", moment->getPopIndices(), factorIds));
              localTriplets.emplace_back(row, col, Scalar(-sign * popIdPower));

              if(popIdPower > 1)
              {
                sslib.dropFactorIds(factorIds, id, 1);

                col = sslib.findCompressedIndex(sslib.getMoment("pi2", moment->getPopIndices(), factorIds));
                localTriplets.emplace_back(row, col, Scalar((popIdPower * (popIdPower - 1)) / 2.));
              }
            }
          }

          else if(countLeft == 1 && countRight == 0)
          {
            if(popIdPower > 0)
            {
              std::vector<size_t> factorIds = moment->getFactorIndices();
              sslib.dropFactorIds(factorIds, id, 1);

              int sign = std::pow(-1, moment->getPopIndices()[0] == id); // sign of contributions that would cancel out if
                                                                        // p1(1-p0) == p0(1-p1)
              col = sslib.findCompressedIndex(sslib.getMoment("pi2", moment->getPopIndices(), factorIds));
              localTriplets.emplace_back(row, col, Scalar(sign * popIdPower));

              if(popIdPower > 1)
              {
                sslib.dropFactorIds(factorIds, id, 1);

                col = sslib.findCompressedIndex(sslib.getMoment("pi2", moment->getPopIndices(), factorIds));
                localTriplets.emplace_back(row, col, Scalar((popIdPower * (popIdPower - 1)) / 2.));
              }
            }
          }

          else if(countLeft == 0 && countRight == 2)
          {
            if(popIdPower > 0)
            {
              std::vector<size_t> factorIds = moment->getFactorIndices();
              sslib.dropFactorIds(factorIds, id, 1);

              col = sslib.findCompressedIndex(sslib.getMoment("Dr", {id, id}, factorIds));
              localTriplets.emplace_back(row, col, Scalar(-popIdPower / 4.));

              factorIds.push_back(sslib.fetchOtherId(id));
              factorIds.push_back(sslib.fetchOtherId(id));

              col = sslib.findCompressedIndex(sslib.getMoment("Dr", {id, id}, factorIds));
              localTriplets.emplace_back(row, col, Scalar(popIdPower / 4.));

              if(popIdPower > 1)
              {
                sslib.dropFactorIds(factorIds, sslib.fetchOtherId(id), 2);
                sslib.dropFactorIds(factorIds, id, 1);

                col = sslib.findCompressedIndex(sslib.getMoment("pi2", moment->getPopIndices(), factorIds));
                localTriplets.emplace_back(row, col, Scalar((popIdPower * (popIdPower - 1)) / 2.));
              }
            }
          }

          else if(countLeft == 0 && countRight == 1)
          {
            if(popIdPower > 0)
            {
              std::vector<size_t> factorIds = moment->getFactorIndices();
              sslib.dropFactorIds(factorIds, id, 1);

              col = sslib.findCompressedIndex(sslib.getMoment("Dr", {id, sslib.fetchOtherId(id)}, factorIds));
              localTriplets.emplace_back(row, col, Scalar(-f * popIdPower / 8.));

              #ifdef NAKED_D
              int signD = std::pow(-1, moment->getPopIndices()[3] == id);
              col = sslib.findCompressedIndex(sslib.getMoment("D", {id}, factorIds));
              localTriplets.emplace_back(row, col, Scalar(-f * signD * (-popIdPower) / 8.));
              #endif

              factorIds.push_back(sslib.fetchOtherId(id));
              factorIds.push_back(sslib.fetchOtherId(id));

              col = sslib.findCompressedIndex(sslib.getMoment("Dr", {id, sslib.fetchOtherId(id)}, factorIds));
              localTriplets.emplace_back(row, col, Scalar(f * popIdPower / 8.));

              #ifdef NAKED_D
              col = sslib.findCompressedIndex(sslib.getMoment("D", {id}, factorIds));
              localTriplets.emplace_back(row, col, Scalar(f * signD * (-popIdPower) / 8.));
              #endif

              if(popIdPower > 1)
              {
                sslib.dropFactorIds(factorIds, id, 1);
                sslib.dropFactorIds(factorIds, sslib.fetchOtherId(id), 2);

                col = sslib.findCompressedIndex(sslib.getMoment("pi2", moment->getPopIndices(), factorIds));
                localTriplets.emplace_back(row, col, Scalar((popIdPower * (popIdPower - 1)) / 2.));
              }
            }
          }

          else if(popIdCount == 0 && popIdPower > 1)
          {
            std::vector<size_t> factorIds = moment->getFactorIndices();
            sslib.dropFactorIds(factorIds, id, 2);

            col = sslib.findCompressedIndex(sslib.getMoment("pi2", moment->getPopIndices(), factorIds));
            localTriplets.emplace_back(row, col, Scalar((popIdPower * (popIdPower - 1)) / 2.));
          }
        } // ends pi2 prefix

        else if(prefix != "I")
          throw bpp::Exception("Drift::mis-specified Moment prefix: " + prefix);
      } // ends parallel loop over basis

      std::vector<Eigen::Triplet<Scalar>> coeffs;
      for(auto& vec : threadTriplets)
        coeffs.insert(coeffs.end(), std::make_move_iterator(vec.begin()), std::make_move_iterator(vec.end()));

      auto matrix = std::make_unique<Matrix<Scalar>>(basisSize, basisSize);
      matrix->setFromTriplets(coeffs);
      matrix->makeCompressed();
      matrix->scale(Scalar(coalRate));

      MatrixEngine::MatrixVariant mv(std::move(*matrix));
      auto engine = std::make_unique<MatrixEngine>(mv, MatrixEngine::VectorVariant{});
      matrices_.emplace_back(std::move(engine));
    } // ends loop over pops
  }
  }, transition_->getMatrixVariant());

  std::cout << "Drift -- assembled individual matrices.\n"
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
      *matrices_[i] *= (newVal / prevVal);
  }

  assembleTransitionMatrix_();
  prevParams_.matchParametersValues(getParameters());
}
