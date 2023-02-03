/* Authors: Gustavo V. Barroso
 * Created: 02/02/2023
 * Last modified: 02/02/2023
 *
 */


#ifndef _DD_MOMENT_H_
#define _DD_MOMENT_H_

#include "Moment.hpp"


class DdMoment:
  public Moment
{

private:
  // ?

public:
  DdMoment():
  Moment()
  { }

  DdMoment(const std::string& name, double value):
  Moment(name, value)
  { }

public:
  // ?

};

#endif
