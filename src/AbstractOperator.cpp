/*
 * Authors: Gustavo V. Barroso
 * Created: 04/04/2023
 * Last modified: 18/09/2025
 *
 */

#include "AbstractOperator.hpp"

void AbstractOperator::printDeltaLDMat(const std::string& fileName)
{
  std::ofstream matFile;
  matFile.open(fileName);

  if(!matFile.is_open())
    throw bpp::Exception("AbstractOperator::failed to open file: " + fileName);

  if(transition_)
  {
    std::visit([&](auto& mat)
    {
      for(int i = 0; i < mat.rows(); ++i)
      {
        for(int j = 0; j < mat.cols(); ++j)
        {
          matFile << mat.eigen().coeff(i, j);

          if(j < mat.cols() - 1)
            matFile << ",";
        }

        matFile << "\n";
      }
    }, transition_->getMatrixVariant());

    matFile.close();
  }

  else
    throw bpp::Exception("AbstractOperator::attempted to print un-initialized transition matrix!");
}

// adds together the different matrices that make up an operator (one per population for Drift;
// population-pair for Migration, etc)
void AbstractOperator::assembleTransitionMatrix_()
{
  MatrixEngine::MatrixVariant combined = matrices_[0]->getMatrixVariant();

  for(size_t i = 1; i < matrices_.size(); ++i)
  {
    combined = std::visit([](auto& a, auto& b) -> MatrixEngine::MatrixVariant
    {
      using AType = std::decay_t<decltype(a)>;
      using BType = std::decay_t<decltype(b)>;

      if constexpr(std::is_same_v<typename AType::Scalar, typename BType::Scalar>)
        return MatrixEngine::MatrixVariant{*a.add(b)};
      else
        throw bpp::Exception("AbstractOperator::add mismatched scalar types between matrices.");
    }, combined, matrices_[i]->getMatrixVariant());
  }


  if(!transition_)
  {
    transition_ = std::make_unique<MatrixEngine>(matrices_[0]->useMPReal());
    transition_->initialize(
        matrices_[0]->getMatrixVariant().index() == 0
            ? std::get<Matrix<double>>(matrices_[0]->getMatrixVariant()).rows()
            : std::get<Matrix<mpfr::mpreal>>(matrices_[0]->getMatrixVariant()).rows(),
        matrices_[0]->getMatrixVariant().index() == 0
            ? std::get<Matrix<double>>(matrices_[0]->getMatrixVariant()).cols()
            : std::get<Matrix<mpfr::mpreal>>(matrices_[0]->getMatrixVariant()).cols(),
        0);
  }

  transition_->setMatrix(combined);
  transition_->addIdentityInPlace(); // convert delta → transition
}
