/*
 * Authors: Gustavo V. Barroso
 * Created: 05/08/2022
 * Last modified: 01/09/2022
 *
 */


#ifndef _SUMSTATSLIBRARY_H_
#define _SUMSTATSLIBRARY_H_

#include <string>
#include <vector>
#include <map>
#include <algorithm>

#include <boost/algorithm/string.hpp>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <Bpp/Text/TextTools.h>

#include "OptionsContainer.hpp"
#include "PolymorphismData.hpp"

class SumStatsLibrary
{

private:
  size_t numPops_;
  size_t order_; // of the moment; number of tracked lineages

  bool initialized_;
  bool compressed_; // indicates whether symmetrical statistics have been pulled together

  // NOTE make this a std::map<std::string, std::pair<size_t, double>> stats_; where the size_t represents the redundancy factor of each statistic (>= 1)?
  // row and column order of matrices follow lexicographical order of stats_'s names
  std::map<std::string, double> stats_;  // name -> value ("observed" Y vector)
  Eigen::Matrix<double, Dynamic, 1> y_; // the Eigen representation

public:
  SumStatsLibrary():
  numPops_(1),
  order_(2),
  initialized_(false),
  compressed_(false),
  stats_(),
  y_()
  { }

  SumStatsLibrary(size_t numPops = 1, size_t order = 2):
  numPops_(numPops),
  order_(order),
  initialized_(false),
  compressed_(false),
  stats_(),
  y_()
  { }

  SumStatsLibrary(const PolymorphismData& dataset):
  numPops_(0),
  order_(0),
  initialized_(false),
  compressed_(false),
  stats_(),
  y_(y)
  {
    init(dataset);
  }

  SumStatsLibrary(const OptionsContainer& opt):
  numPops_(opt.getNumberOfPopulations()),
  order_(opt.getOrder()),
  initialized_(false),
  compressed_(false),
  stats_(),
  y_()
  {
    PolymorphismData dataset;
    dataset.parse(options.getDataPath());
    dataset.computeSumStats();

    init(dataset);
  }

  size_t getNumPops()
  {
    return numPops;
  }


  size_t getOrder()
  {
    return order_;
  }

  bool initialized()
  {
    return initialized_;
  }

  bool compressed()
  {
    return compressed_;
  }

  size_t getNumStats()
  {
    return stats_.size();
  }

  const std::map<std::string, double>& getStats()
  {
    return stats_;
  }

  double getStat(const std::string& name)
  {
    return stats_.at(name);
  }

  const Eigen::Matrix<double, Dynamic, 1>& getYvec()
  {
    return y_;
  }

  void init();

  std::vector<std::string> splitString(const std::string& target, const std::string& query)
  {
    std::vector<std::string> ret(0);
    boost::split(ret, target, boost::is_any_of(query));

    return ret;
  }

  std::string asString(size_t i)
  {
    return bpp::TextTools::toString(i);
  }

  size_t countInstances(const std::string& target, const std::string& query)
  {
    std::string::difference_type count = std::count(std::begin(target), std::end(target), query);
    return static_cast<size_t>(count);
  }

  size_t indexLookup(const std::string& moment)
  {
    auto pos = stats_.find(moment);

    if(pos == std::end(stats_))
      throw bpp::Exception("Moments++::SumStatsLibrary::Could not find index of moment " + moment);

    else
      return static_cast<size_t>(pos - std::begin(stats_));
  }

  void init(size_t numPops, size_t order = 2)
  {
    order_ = order;

    includeHetStats_();
    includeLdStats_();
  }

private:
  void includeHetStats_();

  void includeLdStats_();

  void compress_(); // exploits symmetry among statistics to reduce dimension

};



