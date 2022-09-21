/* Authors: Gustavo V. Barroso
 * Created: 19/09/2022
 * Last modified: 21/09/2022
 *
 */


#ifndef _MOMENT_H_
#define _MOMENT_H_

#include <string>
#include <algorithm>
#include <cstring>
#include <vector>
#include <cstdlib>

#include <boost/algorithm/string.hpp>

#include <Bpp/Text/TextTools.h>


class Moment
{

private:
  std::string name_; // e.g. "Dz_12"
  std::string prefix_; // e.g. "Dz"
  std::vector<size_t> popIndices_; // sorted for binary_search, e.g. {1, 2, 2}
  double value_;

  std::shared_ptr<Moment> parent_; // "equivalent" moment in previous epoch, due to population ancestry

public:
  Moment():
  name_(""),
  prefix_(""),
  popIndices_(0),
  value_(0.),
  parent_(nullptr)
  { }

  Moment(const std::string& name, double value):
  name_(name),
  prefix_(""),
  popIndices_(0),
  value_(value),
  parent_(nullptr)
  {
    std::vector<std::string> splitName(0);
    boost::split(splitName, name, boost::is_any_of("_"));

    prefix_ = splitName[0];

    for(size_t i = 1; i < splitName.size(); ++i)
      popIndices_.push_back(std::stoul(splitName[i]));

    std::sort(std::begin(popIndices_), std::end(popIndices_));
  }

  Moment(const std::string& name, double value, std::shared_ptr<Moment> parent):
  name_(name),
  prefix_(""),
  popIndices_(0),
  value_(value),
  parent_(parent)
  {
    std::vector<std::string> splitName(0);
    boost::split(splitName, name, boost::is_any_of("_"));

    prefix_ = splitName[0];

    for(size_t i = 1; i < splitName.size(); ++i)
      popIndices_.push_back(std::stoul(splitName[i]));

    std::sort(std::begin(popIndices_), std::end(popIndices_));
  }

public:
  void setParent(std::shared_ptr<Moment> mom)
  {
    parent_ = mom;
  }

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

  double getValue() const
  {
    return value_;
  }

  void setValueFromParent()
  {
    value_ = parent_->getValue();
  }

  void setValue(double value)
  {
    value_ = value;
  }

  bool hasPopIndex(size_t index)
  {
    return std::binary_search(std::begin(popIndices_), std::end(popIndices_), index);
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
