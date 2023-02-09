/*
 * Authors: Gustavo V. Barroso
 * Created: 19/09/2022
 * Last modified: 09/02/2023
 *
 */


#ifndef _POPULATION_H_
#define _POPULATION_H_

#include <Bpp/App/ApplicationTools.h>

#include <string>
#include <vector>
#include <memory>
#include <utility>
#include <algorithm>
#include <map>


class Population
{

private:
  std::string name_; // human label, e.g. "Yoruba"
  std::string description_;

  // parents in previous epoch; they will be the same unless focal population is a result of admixture
  std::shared_ptr<Population> leftParent_;
  std::shared_ptr<Population> rightParent_;

  size_t id_; // "i" as it appears in N_i parameters (Drift Operator) and m_ij parameters (Migration Operator)
  size_t startTime_; // in units of generations (coming from past to present)
  size_t endTime_; // in units of generations
  size_t startSize_; // N_i in 1/N_i Drift parameters
  size_t endSize_; // must be equal to startSize because we assume constant pop sizes within each epoch

  bool isDerivedLeftSelected_; // tells q

public:
  Population():
  name_(""),
  description_(""),
  leftParent_(nullptr),
  rightParent_(nullptr),
  id_(0),
  startTime_(0),
  endTime_(0),
  startSize_(0),
  endSize_(0),
  isDerivedLeftSelected_(0)
  { }

  Population(const std::string& name, const std::string& description, size_t id, size_t startTime, size_t endTime, size_t startSize, size_t endSize, bool hasSelection):
  name_(name),
  description_(description),
  leftParent_(nullptr),
  rightParent_(nullptr),
  id_(id),
  startTime_(startTime),
  endTime_(endTime),
  startSize_(startSize),
  endSize_(endSize),
  isDerivedLeftSelected_(hasSelection)
  {
    if(startSize_ != endSize_)
      throw bpp::Exception("Population::start and end sizes for population " + bpp::TextTools::toString(id) + " specified in YAML file are not equal!");
  }

public:
  const std::string& getName() const
  {
    return name_;
  }

  const std::string& getDescription() const
  {
    return description_;
  }

  std::shared_ptr<Population> getLeftParent()
  {
    return leftParent_;
  }

  std::shared_ptr<Population> getRightParent()
  {
    return rightParent_;
  }

  void setLeftParent(std::shared_ptr<Population> parent)
  {
    leftParent_ = parent;
  }

  void setRightParent(std::shared_ptr<Population> parent)
  {
    rightParent_ = parent;
  }

  size_t getId()
  {
    return id_;
  }

  size_t getStartTime()
  {
    return startTime_;
  }

  size_t getEndTime()
  {
    return endTime_;
  }

  bool hasSelection()
  {
    return isDerivedLeftSelected_;
  }

};

#endif
