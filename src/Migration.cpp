/*
 * Authors: Gustavo V. Barroso
 * Created: 10/08/2022
 * Last modified: 10/08/2022
 *
 */


#include "Migration.hpp"


void Migration::setUpMatrix(size_t matrixSize)
{
  matrixSize = 23;
  Eigen::SparseMatrix<int, matrixSize, matrixSize> mat;

  mat(0, 0) = -3;
  mat(1, 0) = 1;
  mat(1, 1) = -1;
  mat(2, 1) = 2;
  mat(2, 2) = -2;
  mat(4, 3) = 1;
  mat(4, 4) = -1;
  mat(5, 3) = 1;
  mat(5, 5) = -1;
  mat(6, 4) = 1;
  mat(6, 5) = 1;
  mat(6, 6) = -2;
  mat(7, 3) = 1;
  mat(7, 7) = -1;
  mat(7, 11) = 4;
  mat(7, 12) = -4;
  mat(7, 14) = -4;
  mat(7, 15) = 4;
  mat(8, 4) = 1;
  mat(8, 7) = 1;
  mat(8, 8) = -2;
  mat(8, 12) = 4;
  mat(8, 13) = -4;
  mat(8, 15) = -4;
  mat(8, 16) = 4;
  mat(9, 5) = 1;
  mat(9, 7) = 1;
  mat(9, 9) = -2;
  mat(9, 14) = 4;
  mat(9, 15) = -4;
  mat(9, 17) = -4;
  mat(9, 18) = 4;
  mat(10, 6) = 1;
  mat(10, 8) = 1;
  mat(10, 9) = 1;
  mat(10, 10) = -3;
  mat(10, 15) = 4;
  mat(10, 16) = -4;
  mat(10, 18) = -4;
  mat(10, 19) = 4;
  mat(12, 11) = 1;
  mat(12, 12) = -1;
  mat(13, 12) = 2;
  mat(13, 13) = -2;
  mat(14, 11) = 1;
  mat(14, 14) = -1;
  mat(15, 12) = 1;
  mat(15, 14) = 1;
  mat(15, 15) = -2;
  mat(16, 13) = 1;
  mat(16, 15) = 2;
  mat(16, 16) = -3;
  mat(17, 14) = 2;
  mat(17, 17) = -2;
  mat(18, 15) = 2;
  mat(18, 17) = 1;
  mat(18, 18) = -3;
  mat(19, 16) = 2;
  mat(19, 18) = 2;
  mat(19, 19) = -4;
  mat(21, 20) = 1;
  mat(21, 21) = -1;
  mat(22, 21) = 2;
  mat(22, 22) = -2;

  matrix_ = mat;

}
