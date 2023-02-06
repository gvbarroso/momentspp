/* Authors: Gustavo V. Barroso
 * Created: 02/02/2023
 * Last modified: 06/02/2023
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
  void printAttributes(std::ostream& stream) override
  {
    stream << std::scientific << position_ << " | " << name_ << " = " << value_ << "; 1st copy derived?: " << isFirstCopyDerived_ << "; constrained?: " << isPutativelySelected_ << ";\n";
  }

  bool isFirstDerived()
  {
    return isFirstCopyDerived_;
  }

  bool isFirstAncestral()
  {
    return !isFirstCopyDerived_;
  }

  bool isSelected()
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
