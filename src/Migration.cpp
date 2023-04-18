/*
 * Authors: Gustavo V. Barroso
 * Created: 10/08/2022
 * Last modified: 31/03/2023
 *
 */


#include "Migration.hpp"

void Migration::setUpMatrices_(const SumStatsLibrary& sslib)
{
  // m_ij is the forward migration rate from pop i to pop j (backwards, the prob that lineage in j has parent in i)
  size_t numPops = littleMigMat_.innerSize();
  size_t sizeOfBasis = sslib.getSizeOfBasis();
  matrices_.reserve(numPops * (numPops - 1));

  for(size_t i = 0; i < numPops; ++i) // for i in m_ij (i => parent pop ID)
  {
    for(size_t j = 0; j < numPops; ++j) // for j in m_ij (j => child pop ID)
    {
      if(i != j)
      {
        std::vector<Eigen::Triplet<double>> coeffs(0);
        coeffs.reserve(3 * sizeOfBasis);

        // although we have pi2(i,j;k,l) moments, we need only to loop over pops twice because we take care of the k,l pop indices in this loop over stats
        for(auto it = std::begin(sslib.getBasis()); it != std::end(sslib.getBasis()); ++it)
        {
          int row = it - std::begin(sslib.getBasis()); // row index
          int col = -1; // inits column index to out-of-bounds
          int childPopIdCount = static_cast<int>((*it)->countInstances(j));

          if(childPopIdCount != 0) // not to populate the sparse matrix with unnecessary zeros
            coeffs.emplace_back(Eigen::Triplet<double>(row, row, -childPopIdCount));

          if((*it)->getPrefix() == "DD")
          {
            std::vector<size_t> popIds = (*it)->getPopIndices();

            for(size_t l = 0; l < popIds.size(); ++ l) // contributions from the DD cols
            {
              if(popIds[l] == j) // if entry matches childPopId
              {
                popIds[l] = i; // assign to focal parentPopId
                col = sslib.findCompressedIndex(sslib.findDdIndex(popIds[0], popIds[1]));
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, 0.5));
                col = sslib.findCompressedIndex(sslib.findDdIndex(popIds[1], popIds[0]));
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, 0.5));
                popIds[l] = j; // recycle
              }
            }

            if((*it)->hasPopIndex(j))
            {
              size_t p1 = popIds[0];
              if(p1 == j)
                p1 = popIds[1];

              size_t p2 = i;
              size_t p3 = j;

              double f = static_cast<double>(p2 == p3) / 2. - 0.25;
              col = sslib.findCompressedIndex(sslib.findDzIndex(p1, p2, p3));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, childPopIdCount * f));
              col = sslib.findCompressedIndex(sslib.findDzIndex(p1, p3, p2));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, childPopIdCount * f));

              p2 = j;

              f = static_cast<double>(p2 == p3) / 2. - 0.25;
              col = sslib.findCompressedIndex(sslib.findDzIndex(p1, p2, p3));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, childPopIdCount * f));

              p2 = i;
              p3 = i;

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
              if(popIds[l] == j) // if entry matches childPopId
              {
                popIds[l] = i; // assign to focal parentPopId
                col = sslib.findCompressedIndex(sslib.findDzIndex(popIds[0], popIds[1], popIds[2]));
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, 0.5));
                col = sslib.findCompressedIndex(sslib.findDzIndex(popIds[0], popIds[2], popIds[1]));
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, 0.5));
                popIds[l] = j; // recycle
              }
            }

            if((*it)->getPopIndices()[0] == j) // contributions from pi2 moments
            {
              // imagine starting with pop indices p2 and p3 on each side of ';' character in pi2(**;**)
              // append parentPopId (i) to the right of both p2 and p3
              // find pi2 statistics by left-right permuting + replacing appendixes

              popIds.clear();
              popIds.resize(4);

              popIds[0] = (*it)->getPopIndices()[1];
              popIds[1] = i;
              popIds[2] = (*it)->getPopIndices()[2];
              popIds[3] = i;

              size_t refCount = std::count(std::begin(popIds), std::end(popIds), j);
              size_t count = std::count(std::begin(popIds), std::end(popIds), j);
              double f = std::pow(-1., count - refCount);

              col = sslib.findCompressedIndex(sslib.findPi2Index(popIds[0], popIds[1], popIds[2], popIds[3]));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, f));
              col = sslib.findCompressedIndex(sslib.findPi2Index(popIds[0], popIds[1], popIds[3], popIds[2]));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, f));
              col = sslib.findCompressedIndex(sslib.findPi2Index(popIds[1], popIds[0], popIds[2], popIds[3]));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, f));
              col = sslib.findCompressedIndex(sslib.findPi2Index(popIds[1], popIds[0], popIds[3], popIds[2]));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, f));

              popIds[3] = j; // switch right appendix to childPopId
              count = std::count(std::begin(popIds), std::end(popIds), j);
              f = std::pow(-1., count - refCount);

              col = sslib.findCompressedIndex(sslib.findPi2Index(popIds[0], popIds[1], popIds[2], popIds[3]));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, f));
              col = sslib.findCompressedIndex(sslib.findPi2Index(popIds[0], popIds[1], popIds[3], popIds[2]));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, f));
              col = sslib.findCompressedIndex(sslib.findPi2Index(popIds[1], popIds[0], popIds[2], popIds[3]));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, f));
              col = sslib.findCompressedIndex(sslib.findPi2Index(popIds[1], popIds[0], popIds[3], popIds[2]));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, f));

              popIds[3] = i; // back
              popIds[1] = j; // switch left appendix to childPopId
              count = std::count(std::begin(popIds), std::end(popIds), j);
              f = std::pow(-1., count - refCount);

              col = sslib.findCompressedIndex(sslib.findPi2Index(popIds[0], popIds[1], popIds[2], popIds[3]));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, f));
              col = sslib.findCompressedIndex(sslib.findPi2Index(popIds[0], popIds[1], popIds[3], popIds[2]));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, f));
              col = sslib.findCompressedIndex(sslib.findPi2Index(popIds[1], popIds[0], popIds[2], popIds[3]));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, f));
              col = sslib.findCompressedIndex(sslib.findPi2Index(popIds[1], popIds[0], popIds[3], popIds[2]));
              coeffs.emplace_back(Eigen::Triplet<double>(row, col, f));

              popIds[3] = j; // have both switched
              count = std::count(std::begin(popIds), std::end(popIds), j);
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

          else if((*it)->getPrefix() == "pi2") // finds pi2 moments that are a single migration away (1-hop neighbors)
          {
            std::vector<size_t> popIds = (*it)->getPopIndices();

            for(size_t l = 0; l < popIds.size(); ++ l) // l -> 0:3
            {
              if(popIds[l] == j) // if entry matches childPopId
              {
                popIds[l] = i; // assign to focal parentPopId
                col = sslib.findCompressedIndex(sslib.findPi2Index(popIds[0], popIds[1], popIds[2], popIds[3]));
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));
                popIds[l] = j; // recycle
              }
            }
          }

          else if((*it)->getPrefix() == "H")
          {
            std::vector<size_t> popIds = (*it)->getPopIndices();

            for(size_t l = 0; l < popIds.size(); ++ l) // l -> 0:1
            {
              if(popIds[l] == j) // if entry matches childPopId
              {
                popIds[l] = i; // assign to focal parentPopId
                col = sslib.findCompressedIndex(sslib.findHetIndex(popIds[0], popIds[1]));
                coeffs.emplace_back(Eigen::Triplet<double>(row, col, 1.));
                popIds[l] = j; // recycle
              }
            }
          }

          else if((*it)->getPrefix() != "I")
            throw bpp::Exception("Migration::mis-specified Moment prefix: " + (*it)->getPrefix());
        }

        Eigen::SparseMatrix<double> mat(sizeOfBasis, sizeOfBasis);
        mat.setFromTriplets(std::begin(coeffs), std::end(coeffs));
        mat.makeCompressed();
        mat *= getParameterValue("m_" + bpp::TextTools::toString(i) + "_" + bpp::TextTools::toString(j));
        matrices_.emplace_back(mat);
      }
    }
  }

  setIdentity_(sizeOfBasis);
  assembleTransitionMatrix_();
}

void Migration::updateMatrices_()
{
  size_t numPops = littleMigMat_.innerSize();
  size_t index = 0;
  std::string paramName = "";

  for(size_t i = 0; i < numPops; ++i)
  {
    for(size_t j = 0; j < numPops; ++j)
    {
      if(i != j)
      {
        paramName = "m_" + bpp::TextTools::toString(i) + "_" + bpp::TextTools::toString(j);

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

void Migration::setLittleMat_()
{
  size_t numPops = fetchNumPops_();
  std::vector<double> row(numPops, 0.);
  Eigen::MatrixXd mat(numPops, numPops);
  mat.setZero();

  for(size_t i = 0; i < numPops; ++i)
  {
    for(size_t j = 0; j < numPops; ++j)
    {
      if(i != j)
        mat(i, j) = getParameterValue("m_" + bpp::TextTools::toString(i) + "_" + bpp::TextTools::toString(j));
    }
  }

  for(size_t i = 0; i < numPops; ++i)
  {
    for(size_t j = 0; j < numPops; ++j)
      row[j] = mat(i, j);

    std::sort(std::begin(row), std::end(row));
    mat(i, i) = 1. - std::accumulate(std::begin(row), std::end(row), 0.);
  }

  littleMigMat_ = mat;
}

void Migration::testFlow_()
{
  int numPops = littleMigMat_.innerSize();
  Graph mig(numPops);

  for(int i = 0; i < numPops; ++i)
  {
    for(int j = 0; j < numPops; ++j)
    {
      if(i != j && littleMigMat_(i, j) > 0.)
        mig.addEdge(i, j);
    }
  }

  for(int i = 0; i < numPops; ++i)
  {
    for(int j = i; j < numPops; ++j)
    {
      bool flow = mig.isReachable(i, j) || mig.isReachable(j, i);

      if(!flow)
      {
        std::cout << "ERROR: Demes " << i << " & " << j << " not connected by migration!\n" << littleMigMat_ << "\n";
        throw bpp::Exception("Migration::no steady-state solution!");
      }
    }
  }
}

size_t Migration::fetchNumPops_() // cute way to get the number of populations P from the raw value of P^2 - P ( == matrices_.size())
{
  int numPops = 2; // we want the positive solution of the quadratic equation P^2 - P - matrices_.size() = 0
  int n = static_cast<int>(getParameters().size()); // raw value of P^2 - P

  for(int i = 2; i < n; ++i)
  {
    if(i * (1 - i) == -n)  // guaranteed to find if matrices_.size() was built correctly
    {
      numPops = i;
      break;
    }
  }

  return static_cast<size_t>(numPops);
}

