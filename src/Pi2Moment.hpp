/* Authors: Gustavo V. Barroso
 * Created: 02/02/2023
 * Last modified: 06/02/2023
 *
 */


#ifndef _PI2_MOMENT_H_
#define _PI2_MOMENT_H_

#include <memory>

#include "Moment.hpp"
#include "HetMoment.hpp"

class Pi2Moment: public Moment
{

private:
  std::shared_ptr<HetMoment> left_; // left locus
  std::shared_ptr<HetMoment> right_; // right locus

public:
  Pi2Moment():
  Moment(),
  left_(nullptr),
  right_(nullptr)
  { }

  Pi2Moment(const std::string& name, double value, std::shared_ptr<HetMoment> left, std::shared_ptr<HetMoment> right):
  Moment(name, value),
  left_(left),
  right_(right)
  { }

public:
  void printAttributes(std::ostream& stream) override
  {
    stream << std::scientific << position_ << " | " << name_ << " = " << value_ << "\nleftHet: ";
    left_->printAttributes(stream);
    stream << "rightHet: ";
    right_->printAttributes(stream);
    stream << "\n\n";
  }

  std::shared_ptr<HetMoment> getLeftHetStat()
  {
    return left_;
  }

  std::shared_ptr<HetMoment> getRightHetStat()
  {
    return right_;
  }

  void setLeftHetStat(std::shared_ptr<HetMoment> mom)
  {
    left_ = mom;
  }

  void setRightHetStat(std::shared_ptr<HetMoment> mom)
  {
    right_ = mom;
  }

};

#endif
