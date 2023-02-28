/* Authors: Gustavo V. Barroso
 * Created: 19/09/2022
 * Last modified: 13/02/2023
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

protected:
  std::string name_; // e.g. "Dz_12_A"
  std::string prefix_; // e.g. "Dz"
  std::string suffix_; // e.g. "A", refers to the permutation of sampling order of derived/ancestral (H and Pi2 stats)
  std::vector<size_t> popIndices_; // sorted for binary_search, e.g. {1, 2, 2}, may be relevant if Order and numPops are high
  size_t position_; // within the Y vector, in SumStatsLibrary, same as row/col in AbstractOperator matrix(ces)
  double value_;

  std::shared_ptr<Moment> parent_; // "equivalent" moment in previous epoch, according to population ancestry (via popIndices_)
  std::vector<std::shared_ptr<Moment>> aliases_; // equivalent moments (permuations with same expectations) under a given scenario (eg selection against derived allele in left locus)

public:
  Moment():
  name_(""),
  prefix_(""),
  suffix_(""),
  popIndices_(0),
  position_(0),
  value_(0.),
  parent_(nullptr),
  aliases_(0)
  { }

  Moment(const std::string& name, double value):
  name_(name),
  prefix_(""),
  suffix_(""),
  popIndices_(0),
  position_(0),
  value_(value),
  parent_(nullptr),
  aliases_(0)
  {
    parseName_(name);
  }

public:
  virtual ~Moment() = default;

  virtual void printAttributes(std::ostream& stream)
  {
    stream << std::scientific << position_ << " | " << name_ << " = " << value_;

    if(parent_ != nullptr)
      stream << "; parent = " << parent_->getName() << "\n";

    else
      stream << "\n";

    if(aliases_.size() > 0)
    {
      stream << "\taliases: ";

      for(auto it = std::begin(aliases_); it != std::end(aliases_); ++it)
        stream << (*it)->getName() << ",";

      stream << "\n";
    }
  }

  const std::string& getName() const
  {
    return name_;
  }

  const std::string& getPrefix() const
  {
    return prefix_;
  }

  const std::string& getSuffix() const
  {
    return suffix_;
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

  const std::vector<std::shared_ptr<Moment>>& getAliases()
  {
    return aliases_;
  }

  size_t getNumberOfAliases()
  {
    return aliases_.size();
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

  void insertAlias(std::shared_ptr<Moment> mom)
  {
    if(std::find(std::begin(aliases_), std::end(aliases_), mom) == std::end(aliases_))
      aliases_.push_back(mom);

    else
      throw bpp::Exception("Moment::attempted to duplicate alias " + name_ + "->" + mom->getName());
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

  bool hasAlias(std::shared_ptr<Moment> mom)
  {
    return std::find(std::begin(aliases_), std::end(aliases_), mom) != std::end(aliases_);
  }

private:
  void parseName_(const std::string& name)
  {
    std::vector<std::string> splitName(0);
    boost::split(splitName, name, boost::is_any_of("_"));
    prefix_ = splitName[0];

    if(splitName.size() > 1)
    {
      suffix_ = splitName.back();

      for(size_t i = 1; i < splitName.size() - 1; ++i)
        popIndices_.emplace_back(std::stoul(splitName[i]));
    }
  }
};

#endif
