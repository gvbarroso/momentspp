/* Authors: Gustavo V. Barroso
 * Created: 19/09/2022
 * Last modified: 03/06/2024
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
#include <typeinfo>

#include <boost/algorithm/string.hpp>

#include <Bpp/Text/TextTools.h>
#include <Bpp/Exceptions.h>


class Moment
{

protected:
  std::string name_; // e.g. "pi2_1100_(1-2p0)^2"; within each Epoch, name is unique to *this Moment
  std::string prefix_; // e.g. "pi2"
  std::vector<size_t> popIndices_; // population indices associated with main statistic
  std::vector<size_t> factorIndices_; // population indices associated with each (1-2p) factor
  size_t position_; // index within the Y vector (see Epoch::computeExpectedSumStats()) and SumStatsLibrary basis_
  long double value_; // expectation

  std::shared_ptr<Moment> parent_; // "equivalent" moment in previous epoch, according to population ancestry
  std::vector<std::weak_ptr<Moment>> aliases_; // equivalent moments (permutations with same expectations)

public:
  Moment():
  name_(""),
  prefix_(""),
  popIndices_(0),
  factorIndices_(0),
  position_(0),
  value_(0.),
  parent_(nullptr),
  aliases_(0)
  { }

  Moment(const std::string& name, long double value):
  name_(name),
  prefix_(""),
  popIndices_(0),
  factorIndices_(0),
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
  }

  virtual bool hasSamePopIds(std::shared_ptr<Moment> mom)
  {
    return mom->getPopIndices() == popIndices_;
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

  const std::vector<size_t>& getFactorIndices() const
  {
    return factorIndices_;
  }

  size_t getFactorPower()
  {
    return factorIndices_.size();
  }

  int getPopFactorPower(size_t index)
  {
    return std::count(std::begin(factorIndices_), std::end(factorIndices_), index);
  }

  int getPopFactorPower(size_t index) const
  {
    return std::count(std::begin(factorIndices_), std::end(factorIndices_), index);
  }

  size_t getPosition() const
  {
    return position_;
  }

  long double getValue() const
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
    if(popIndices_.size() == 1) // naked signed D
      return false;

    else
    {
      auto tmp = popIndices_;
      std::sort(std::begin(tmp), std::end(tmp));

      return std::adjacent_find(std::begin(popIndices_), std::end(popIndices_), std::not_equal_to<size_t>()) != std::end(popIndices_);
    }
  }

  bool hasCrossPopFactors()
  {
    if(factorIndices_.size() > 0)
      return std::adjacent_find(std::begin(factorIndices_), std::end(factorIndices_), std::not_equal_to<size_t>()) != std::end(factorIndices_);

    else
      return false;
  }

  void setValueFromParent()
  {
    value_ = parent_->getValue();
  }

  void setValue(long double value)
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

  bool hasAnyOfPopIndices(const std::vector<size_t>& ids) const
  {
    for(size_t i = 0; i < ids.size(); ++i)
    {
      if(!(std::find(std::begin(popIndices_), std::end(popIndices_), ids[i]) == std::end(popIndices_)))
        return true;
    }

    return false;
  }

  bool hasAllOfPopIndices(const std::vector<size_t>& ids) const
  {
    for(size_t i = 0; i < ids.size(); ++i)
    {
      if((std::find(std::begin(popIndices_), std::end(popIndices_), ids[i]) == std::end(popIndices_)))
        return false;
    }

    return true;
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
                            { return tmp.lock()->getName() == obj.lock()->getName(); });

    return pos != std::end(aliases_);
  }

  std::vector<size_t> fetchDiffPopIds(size_t focalPopId) // fetches pop IDs different from focalPopId
  {
    std::vector<size_t> ret(0);

    for(auto it = std::begin(popIndices_); it != std::end(popIndices_); ++it)
    {
      if(*it != focalPopId)
        ret.push_back(*it);
    }

    return ret;
  }

  size_t popIndexDistance(const std::shared_ptr<Moment> other)
  {
    static_assert(std::is_same_v<decltype(*this), decltype(*other.get())>);

    size_t dist = 0;
    for(size_t i = 0; i < popIndices_.size(); ++i)
    {
      if(popIndices_[i] != other->getPopIndices()[i])
        ++dist;
    }

    return dist;
  }

  bool isMigAdjacent(const std::shared_ptr<Moment> other) // NOTE check
  {
    return popIndexDistance(other) <= 1;
  }

  // directional, tells if *this can be reached by other* via Admixture
  bool isAdmixAdjacent(const std::shared_ptr<Moment> other, size_t fromId, size_t toId)
  {
    if(!std::is_same_v<decltype(*this), decltype(*other.get())>)
      return 0;

    bool eval = 1;

    for(size_t i = 0; i < popIndices_.size(); ++i)
    {
      if(popIndices_[i] != other->getPopIndices()[i])
      {
        if(popIndices_[i] != toId)
          eval = 0;

        else if(other->getPopIndices()[i] != fromId)
          eval = 0;
      }
    }

    return eval;
  }

  void parseName(const std::string& name)
  {
    std::vector<std::string> splitName(0);
    boost::split(splitName, name, boost::is_any_of("_"));
    prefix_ = splitName[0];

    if(splitName.size() > 1)
    {
      for(size_t i = 1; i < splitName.size(); ++i)
      {
        if(splitName[i] != "l")
          popIndices_.push_back(std::stoul(splitName[i]));

        else
        {
          for(size_t j = i + 1; j < splitName.size(); ++j)
            factorIndices_.push_back(std::stoul(splitName[j]));

          break;
        }
      }
    }

    // sorts to compare *this and other moments more easily NOTE popIndices_ should not be sorted!
    std::sort(std::begin(factorIndices_), std::end(factorIndices_));

    // aesthetics
    std::string nome = prefix_;
    for(size_t i = 0; i < popIndices_.size(); ++i)
      nome = nome + "_" + bpp::TextTools::toString(popIndices_[i]);

    std::vector<size_t> indices = factorIndices_;
    indices.erase(std::unique(std::begin(indices), std::end(indices)), std::end(indices));

    for(size_t i = 0; i < indices.size(); ++i)
    {
      size_t count = std::count(std::begin(factorIndices_), std::end(factorIndices_), indices[i]);

      if(count > 0)
        nome = nome + "_(1-2p" + bpp::TextTools::toString(indices[i]) + ")^" + bpp::TextTools::toString(count);
    }

    name_ = nome;
  }
};

#endif
