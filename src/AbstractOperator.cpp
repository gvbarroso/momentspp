/*
 * Authors: Gustavo V. Barroso
 * Created: 09/08/2022
 * Last modified: 24/09/2025
 *
 */


#include "AbstractOperator.hpp"

#include <fstream>

void AbstractOperator::printDeltaLDMat(const std::string& fileName)
{
  std::ofstream matFile(fileName);

  if(!matFile.is_open())
    throw bpp::Exception("AbstractOperator::failed to open file: " + fileName);

  if(!transition_)
    throw bpp::Exception("AbstractOperator::attempted to print un-initialized transition matrix!");

  auto eigenVar = transition_->toEigenMatrixVariant();

  std::visit(overloaded{[&](const auto& mat)
  {
    const int rows = mat.rows();
    const int cols = mat.cols();

    for(int i = 0; i < rows; ++i)
    {
      for(int j = 0; j < cols; ++j)
      {
        matFile << mat.coeff(i, j);

        if(j + 1 < cols)
          matFile << ",";
      }

      matFile << "\n";
    }
  }
  }, eigenVar);

  matFile.close();
}
