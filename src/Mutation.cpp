/*
 * Authors: Gustavo V. Barroso
 * Created: 10/08/2022
 * Last modified: 10/08/2022
 *
 */


#include "Mutation.hpp"


void Mutation::setUpMatrix(size_t matrixSize)
{
  matrixSize = 23;
  Eigen::SparseMatrix<int, matrixSize, matrixSize> mat;

  mat(0, 0) = -3;
  mat(0, 4) = 1;
  mat(1, 1) = -1;
  mat(3, 0) = 4;
  mat(3, 3) = -5;
  mat(4, 4) = -3;
  mat(5, 5) = -3;
  mat(6, 6) = -1;
  mat(11, 3) = 1;
  mat(11, 11) = -2;
  mat(12, 12) = -1;
  mat(13, 13) = -1;
  mat(14, 14) = -1;
  mat(18, 18) = -1;
  mat(21, 21) = -3;

  matrix_ = mat;

}
