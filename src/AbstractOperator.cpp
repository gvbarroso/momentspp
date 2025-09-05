/*
 * Authors: Gustavo V. Barroso
 * Created: 04/04/2023
 * Last modified: 05/09/2025
 *
 */


#include "AbstractOperator.hpp"


void AbstractOperator::printDeltaLDMat(const std::string& fileName)
{
  std::ofstream matFile;
  matFile.open(fileName);

  auto mat = matrices_[0];

  if(matrices_.size() > 1)
  {
    for(size_t i = 1; i < matrices_.size(); ++i)
      mat += matrices_[i];
  }

  for(int i = 0; i < mat.rows(); ++i)
  {
    for(int j = 0; j < mat.cols(); ++j)
    {
      matFile << mat.coeffRef(i, j);

      if(j < mat.cols() - 1)
        matFile << ",";
    }

    matFile  << "\n";
  }

  matFile.close();
}

void AbstractOperator::printTransitionLDMat(const std::string& fileName)
{
  std::ofstream matFile;
  matFile.open(fileName);

  for(int i = 0; i < transition_.rows(); ++i)
  {
    for(int j = 0; j < transition_.cols(); ++j)
    {
      matFile << transition_.coeffRef(i, j);

      if(j < transition_.cols() - 1)
        matFile << ",";

    }

    matFile  << "\n";
  }

  matFile.close();
}

// adds together the different matrices that make up an operator (one per population for Drift; population-pair for Migration, etc)
void AbstractOperator::assembleTransitionMatrix_()
{
  transition_ = matrices_[0]; // inits to "delta" matrix

  if(matrices_.size() > 1)
  {
    for(size_t i = 1; i < matrices_.size(); ++i)
      transition_ += matrices_[i];
  }
  
  // converts from "delta" to "transition" matrix
  //transition_ += identity_; // NOTE: not used ATM because we are SUMMING matrices from each Operator in Epoch::init_()
}

void AbstractOperator::setIdentity_(size_t numStats)
{
  Eigen::SparseMatrix<mpfr::mpreal> id(numStats, numStats);

  std::vector<Eigen::Triplet<mpfr::mpreal>> md(0);
  md.reserve(numStats);

  for(size_t i = 0; i < numStats; ++i)
    md.emplace_back(Eigen::Triplet<mpfr::mpreal>(i, i, 1.));

  id.setFromTriplets(std::begin(md), std::end(md));
  id.makeCompressed();

  identity_ = id;
}
