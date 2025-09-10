/*
 * Authors: Gustavo V. Barroso
 * Created: 04/04/2023
 * Last modified: 09/09/2025
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
  // clones / inits to "delta" matrix
  std::unique_ptr<MatrixInterface> sum = matrices_[0]->add(*matrices_[0]); // Identity operation

  if(matrices_.size() > 1)
  {
    for(size_t i = 1; i < matrices_.size(); ++i)
      sum = sum->add(*matrices_[i]);
  }

  transition_ = sum;
}

void AbstractOperator::setIdentity_(size_t numStats)
{
  std::vector<Eigen::Triplet<double>> md;
  md.reserve(numStats);

  for(size_t i = 0; i < numStats; ++i)
    md.emplace_back(i, i, 1.0);

  if(!identity_)
  {
    if(dynamic_cast<MatrixDouble*>(matrices_[0].get()))
      identity_ = std::make_unique<MatrixDouble>(numStats, numStats);

    else if(dynamic_cast<MatrixMPReal*>(matrices_[0].get()))
      identity_ = std::make_unique<MatrixMPReal>(numStats, numStats);

    else
      throw bpp::Exception("AbstractOperator::Mis-cast transition matrix!");
  }

  identity_->setFromTriplets(md);
  identity_->makeCompressed();
}
