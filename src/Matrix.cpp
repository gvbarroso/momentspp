/*
 * Authors: Gustavo V. Barroso
 * Created: 28/07/2022
 * Last modified: 28/02/2022
 *
 */

/*
 * This file is modified after code from the zipHMM libraries:
 * https://github.com/mailund/ziphmm
 */



#include "matrix.hpp"

std::ostream &operator<<(std::ostream &out, const Matrix& mat) {
  for(unsigned i = 0; i < mat.get_height(); ++i)
  {
    for(unsigned j = 0; j < mat.get_width(); ++j)
      out << mat(i, j) << " ";
    out << std::endl;
  }
  
  return out;
}

