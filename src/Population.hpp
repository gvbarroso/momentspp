/*
 * Authors: Gustavo V. Barroso
 * Created: 19/09/2022
 * Last modified: 20/09/2022
 *
 */


#ifndef _POPULATION_H_
#define _POPULATION_H_

#include <string>
#include <vector>
#include <memory>
#include <utility>
#include <algorithm>
#include <map>


class Population
{
private:
  size_t id_; // "i" as it appears in N_i parameters (Drift Operator) and m_ij parameters (Migration Operator)
  std::string name_; // human label, e.g. "Yoruba"

  // parents in previous epoch; they will be the same unless focal population is a result of admixture
  std::shared_ptr<Population> leftParent_;
  std::shared_ptr<Population> rightParent_;

public:
  Population(size_t id, const std::string& name):
  id_(id),
  name_(name),
  leftParent_(nullptr),
  rightParent_(nullptr),
  stats_()
  { }

public:
  size_t getId()
  {
    return id_;
  }

  const std::string& getName() const
  {
    return name_;
  }

  std::shared_ptr<Population> getLeftParent()
  {
    return leftParent_;
  }

  std::shared_ptr<Population> getRightParent()
  {
    return rightParent_;
  }

};

#endif
