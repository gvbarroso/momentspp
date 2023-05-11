/* Authors: Gustavo V. Barroso
 * Created: 02/02/2023
 * Last modified: 03/03/2023
 *
 */


#ifndef _HET_MOMENT_H_
#define _HET_MOMENT_H_

#include "Moment.hpp"


class HetMoment: public Moment
{

private:
  bool isPutativelySelected_;

public:
  HetMoment():
  Moment(),
  isPutativelySelected_(0)
  { }

  HetMoment(const std::string& name, double value, bool isPutativelySelected):
  Moment(name, value),
  isPutativelySelected_(isPutativelySelected)
  { }

public:
  bool isConstrained()
  {
    return isPutativelySelected_;
  }

  bool isNeutral()
  {
    return !isPutativelySelected_;
  }

  void setConstraintStatus(bool isSelected)
  {
    isPutativelySelected_ = isSelected;
  }

  bool isCrossPop()
  {
    return popIndices_[0] != popIndices_[1];
  }

  bool hasSamePopIds(std::shared_ptr<Moment> other) override
  {
    assert(std::dynamic_pointer_cast<HetMoment>(other) != nullptr);

    bool test = 0;

    if(other->getPopIndices() == popIndices_)
      test = 1;

    else if(other->getPopIndices()[1] == popIndices_[0] && other->getPopIndices()[0] == popIndices_[1])
      test = 1;

    return test;
  }

};

#endif
