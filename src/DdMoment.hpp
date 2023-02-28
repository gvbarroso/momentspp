/* Authors: Gustavo V. Barroso
 * Created: 02/02/2023
 * Last modified: 27/02/2023
 *
 */


#ifndef _DD_MOMENT_H_
#define _DD_MOMENT_H_

#include "Moment.hpp"


class DdMoment: public Moment
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

  bool hasSamePopIds(std::shared_ptr<Moment> mom) override
  {
    assert(std::dynamic_pointer_cast<DdMoment>(mom) != nullptr);

    bool test = 0;

    if(mom->getPopIndices() == popIndices_)
      test = 1;

    else if(mom->getPopIndices()[1] == popIndices_[0] && mom->getPopIndices()[0] == popIndices_[1])
      test = 1;

    return test;
  }

};

#endif
