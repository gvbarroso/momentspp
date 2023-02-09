/* Authors: Gustavo V. Barroso
 * Created: 02/02/2023
 * Last modified: 07/02/2023
 *
 */


#ifndef _HET_MOMENT_H_
#define _HET_MOMENT_H_

#include "Moment.hpp"


class HetMoment: public Moment
{

private:
  bool isFirstCopyDerived_;
  bool isPutativelySelected_;

public:
  HetMoment():
  Moment(),
  isFirstCopyDerived_(0),
  isPutativelySelected_(0)
  { }

  HetMoment(const std::string& name, double value, bool isFirstDerived, bool isPutativelySelected):
  Moment(name, value),
  isFirstCopyDerived_(isFirstDerived),
  isPutativelySelected_(isPutativelySelected)
  { }

public:
  bool isFirstDerived()
  {
    return isFirstCopyDerived_;
  }

  bool isFirstAncestral()
  {
    return !isFirstCopyDerived_;
  }

  bool isConstrained()
  {
    return isPutativelySelected_;
  }

  bool isNeutral()
  {
    return !isPutativelySelected_;
  }

  void setFirstCopyStatus(bool isDerived)
  {
    isFirstCopyDerived_ = isDerived;
  }

  void setConstraintStatus(bool isSelected)
  {
    isPutativelySelected_ = isSelected;
  }

};

#endif
