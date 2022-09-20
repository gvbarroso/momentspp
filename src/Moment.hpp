/* Authors: Gustavo V. Barroso
 * Created: 19/09/2022
 * Last modified: 19/09/2022
 *
 */


#ifndef _MOMENT_H_
#define _MOMENT_H_

#include <string>
#include <algorithm>
#include <cstring>
#include <vector>

#include <boost/algorithm/string.hpp>

#include <Bpp/Text/TextTools.h>


class Moment
{

private:
  std::string name_; // e.g. "Dz_12"
  std::string prefix_; // e.g. "Dz"
  std::vector<size_t> popIndices_; // sorted, e.g. {1, 2}
  double value_;

public:
  Moment(const std::string& name, double value):
  name_(name),
  popIndices_(0),
  value_(value)
  {

    std::sort(std::begin(popIndices_), std::end(popIndices_));
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

  double getValue()
  {
    return value_;
  }

  void setValue(double value)
  {
    value_ = value;
  }

  bool hasPopIndex(size_t index)
  {
    return std::binary_search(std::begin(popIndices_), std::end(popIndices_), index);
  }

};

#endif
