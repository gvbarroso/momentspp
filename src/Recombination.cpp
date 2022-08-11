/*
 * Authors: Gustavo V. Barroso
 * Created: 09/08/2022
 * Last modified: 09/08/2022
 *
 */


#include "Recombination.hpp"


void Recombination::setUpMatrix(size_t matrixSize)
{
  matrixSize = 23;
  Eigen::SparseMatrix<int, matrixSize, matrixSize> mat;

  mat(0, 0) = -2;
  mat(1, 1) = -2;
  mat(2, 2) = -2;
  mat(3, 3) = -1;
  mat(4, 4) = -1;
  mat(5, 5) = -1;
  mat(6, 6) = -1;
  mat(7, 7) = -1;
  mat(8, 8) = -1;
  mat(9, 9) = -1;
  mat(10, 10) = -1;

  matrix_ = mat;

}
