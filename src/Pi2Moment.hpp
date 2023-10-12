/* Authors: Gustavo V. Barroso
 * Created: 02/02/2023
 * Last modified: 12/10/2023
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

  void printAttributes(std::ostream& stream) override
  {
    stream << "(" << position_ << ") " << name_;
    stream << " = " << std::scientific << value_;

    if(parent_ != nullptr)
      stream << "; parent = " << parent_->getName() << "\n";

    else
      stream << "\n";

    if(aliases_.size() > 0)
    {
      auto tmp = getAliases();
      stream << "\talias(es): ";

      for(auto it = std::begin(tmp); it != std::end(tmp); ++it)
      {
        stream << (*it)->getName();
        auto test = it;

        if(++test != std::end(tmp))
          stream << ",";
      }

      stream << "\n";
    }

    stream << "L = " << left_ << " ";
    left_->printAttributes(stream);

    stream << "R = " << right_ << " ";
    right_->printAttributes(stream);
  }

  bool hasSamePopIds(std::shared_ptr<Moment> other) override
  {
    assert(std::dynamic_pointer_cast<Pi2Moment>(other) != nullptr);

    bool test = 0;

    if(other->getPopIndices() == popIndices_)
      test = 1;

    else
    {
      auto tmpThis = popIndices_;
      auto tmpOther = other->getPopIndices();

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
