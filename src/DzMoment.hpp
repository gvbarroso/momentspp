/* Authors: Gustavo V. Barroso
 * Created: 02/02/2023
 * Last modified: 02/02/2023
 *
 */


#ifndef _DZ_MOMENT_H_
#define _DZ_MOMENT_H_

#include "Moment.hpp"


class DzMoment:
  public Moment
{

private:
  // ?

public:
  DzMoment():
  Moment()
  { }

  DzMoment(const std::string& name, double value):
  Moment(name, value)
  { }

public:
  // ?

};

#endif
