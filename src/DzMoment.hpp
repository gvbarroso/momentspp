/* Authors: Gustavo V. Barroso
 * Created: 02/02/2023
 * Last modified: 27/02/2023
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

  bool hasSamePopIds(std::shared_ptr<Moment> mom) override
  {
    assert(std::dynamic_pointer_cast<DzMoment>(mom) != nullptr);

    bool test = 0;

    if(mom->getPopIndices() == popIndices_)
      test = 1;

    else
    {
      auto tmpThis = popIndices_;
      auto tmpOther = mom->getPopIndices();

      std::sort(std::begin(tmpThis), std::end(tmpThis));
      tmpThis.erase(std::unique(std::begin(tmpThis), std::end(tmpThis)), std::end(tmpThis));

      std::sort(std::begin(tmpOther), std::end(tmpOther));
      tmpOther.erase(std::unique(std::begin(tmpOther), std::end(tmpOther)), std::end(tmpOther));

      if(tmpThis == tmpOther)
        test = 1;
    }

    return test;
  }

};

#endif
