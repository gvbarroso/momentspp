/* Authors: Gustavo V. Barroso
 * Created: 19/09/2022
 * Last modified: 03/04/2023
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
#include <cassert>
#include <memory>

#include <boost/algorithm/string.hpp>

#include <Bpp/Text/TextTools.h>


class Moment
{

protected:
  std::string name_; // e.g. "Dz_12_A"
  std::string prefix_; // e.g. "Dz"
  std::vector<size_t> popIndices_; // sorted for binary_search, e.g. {1, 2, 2}, may be relevant if Order and numPops are high
  size_t position_; // within the Y vector, in SumStatsLibrary, same as row/col in AbstractOperator matrix(ces)
  double value_;

  std::shared_ptr<Moment> parent_; // "equivalent" moment in previous epoch, according to population ancestry (via popIndices_)
  std::vector<std::weak_ptr<Moment>> aliases_; // equivalent moments (permuations with same expectations) under a given scenario (eg selection against derived allele in left locus)

public:
  Moment():
  name_(""),
  prefix_(""),
  popIndices_(0),
  position_(0),
  value_(0.),
  parent_(nullptr),
  aliases_(0)
  { }

  Moment(const std::string& name, double value):
  name_(name),
  prefix_(""),
  popIndices_(0),
  position_(0),
  value_(value),
  parent_(nullptr),
  aliases_(0)
  {
    parseName(name);
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
      auto tmp = getAliases();
      stream << "\taliases: ";

      for(auto it = std::begin(tmp); it != std::end(tmp); ++it)
        stream << (*it)->getName() << ",";

      stream << "\n";
    }
  }

  virtual bool hasSamePopIds(std::shared_ptr<Moment> mom)
  {
    return mom->getPopIndices() == popIndices_; // default, needed for DummyMoment "I" (Base)
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

  std::vector<std::shared_ptr<Moment>> getAliases()
  {
    std::vector<std::shared_ptr<Moment>> tmp(0);
    tmp.reserve(aliases_.size());

    for(auto it = std::begin(aliases_); it != std::end(aliases_); ++it)
      tmp.emplace_back((*it).lock());

    return tmp;
  }

  size_t getNumberOfAliases()
  {
    return aliases_.size();
  }

  bool isCrossPop() // tells if moment involves samples from more than 1 population
  {
    return std::adjacent_find(std::begin(popIndices_), std::end(popIndices_), std::not_equal_to<size_t>()) != std::end(popIndices_);
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
    if(!hasAlias(mom))
    {
      std::weak_ptr<Moment> tmp = mom;
      aliases_.push_back(tmp);
    }

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
    std::weak_ptr<Moment> tmp = mom;

    auto pos = std::find_if(std::begin(aliases_), std::end(aliases_), [&tmp](const auto& obj)
                            { return tmp.lock() == obj.lock(); });

    return pos != std::end(aliases_);
  }

  std::vector<size_t> fetchDiffPopIds(size_t notPopId) // fetches pop IDs different from notPopId
  {
    std::vector<size_t> ret(0);

    for(auto it = std::begin(popIndices_); it != std::end(popIndices_); ++it)
    {
      if(*it != notPopId)
        ret.push_back(*it);
    }

    return ret;
  }

  void parseName(const std::string& name)
  {
    std::vector<std::string> splitName(0);
    boost::split(splitName, name, boost::is_any_of("_"));
    prefix_ = splitName[0];

    if(splitName.size() > 1)
    {
      for(size_t i = 1; i < splitName.size(); ++i)
        popIndices_.push_back(std::stoul(splitName[i]));
    }

    //std::sort(std::begin(popIndices_), std::end(popIndices_));
  }

};

#endif
