/* Authors: Gustavo V. Barroso
 * Created: 02/02/2023
 * Last modified: 09/02/2023
 *
 */


#ifndef _DZ_MOMENT_H_
#define _DZ_MOMENT_H_

#include "Moment.hpp"


class DzMoment: public Moment
{

private:
  bool isConstrained_; // does j in Dz_i_j_k concern a population with selection on the derived allele?

public:
  DzMoment():
  Moment(),
  isConstrained_(0)
  { }

  DzMoment(const std::string& name, double value):
  Moment(name, value),
  isConstrained_(0)
  { }

public:
  bool isConstrained()
  {
    return isConstrained_;
  }

};

#endif
