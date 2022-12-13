/* Authors: Gustavo V. Barroso
 * Created: 19/09/2022
 * Last modified: 13/12/2022
 *
 */


#ifndef _MOMENT_H_
#define _MOMENT_H_

#include <string>
#include <algorithm>
#include <cstring>
#include <vector>
#include <cstdlib>
#include <iostream>

#include <boost/algorithm/string.hpp>

#include <Bpp/Text/TextTools.h>


class Moment
{

private:
  std::string name_; // e.g. "Dz_12"
  std::string prefix_; // e.g. "Dz"
  std::vector<size_t> popIndices_; // sorted for binary_search, e.g. {1, 2, 2}
  size_t position_; // in the Y vector
  double value_;

  std::shared_ptr<Moment> parent_; // "equivalent" moment in previous epoch, according to population ancestry

public:
  Moment():
  name_(""),
  prefix_(""),
  popIndices_(0),
  position_(0),
  value_(0.),
  parent_(nullptr)
  { }

  Moment(const std::string& name, double value):
  name_(name),
  prefix_(""),
  popIndices_(0),
  position_(0),
  value_(value),
  parent_(nullptr)
  {
    std::vector<std::string> splitName(0);
    boost::split(splitName, name, boost::is_any_of("_"));
    prefix_ = splitName[0];

    if(splitName.size() > 1)
    {
      for(size_t i = 1; i < splitName.size(); ++i)
        popIndices_.push_back(std::stoul(splitName[i]));
    }
  }

public:
  const std::string& getName() const
  {
    return name_;
  }

  const std::string& getPrefix() const
  {
    return prefix_;
  }

  const std::vector<size_t>& getPopIndices() const
  {
    return popIndices_;
  }

  size_t getPosition() const
  {
    return position_;
  }

  double getValue() const
  {
    return value_;
  }

  const std::shared_ptr<Moment> getParent()
  {
    return parent_;
  }

  void setValueFromParent()
  {
    value_ = parent_->getValue();
  }

  void setValue(double value)
  {
    value_ = value;
  }

  void setParent(std::shared_ptr<Moment> mom)
  {
    parent_ = mom;
  }

  void setPosition(size_t pos)
  {
    position_ = pos;
  }

  bool hasPopIndex(size_t index)
  {
    return !(std::find(std::begin(popIndices_), std::end(popIndices_), index) == std::end(popIndices_));
  }

  size_t countInstances(size_t index) const
  {
    return std::count(std::begin(popIndices_), std::end(popIndices_), index);
  }

  size_t countInstances(size_t index)
  {
    return std::count(std::begin(popIndices_), std::end(popIndices_), index);
  }

};

#endif
