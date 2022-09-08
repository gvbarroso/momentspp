/*
 * Authors: Gustavo V. Barroso
 * Created: 05/08/2022
 * Last modified: 08/09/2022
 *
 */


#ifndef _SUMSTATSLIBRARY_H_
#define _SUMSTATSLIBRARY_H_

#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <cstring>

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
  Eigen::VectorXd y_; // the Eigen representation of the observed vector
  Eigen::MatrixXd covar_; // covariance matrix of observed sum stats

public:
  SumStatsLibrary():
  numPops_(1),
  order_(2),
  initialized_(false),
  compressed_(false),
  stats_(),
  y_(),
  covar_()
  { }

  SumStatsLibrary(size_t numPops = 1, size_t order = 2):
  numPops_(numPops),
  order_(order),
  initialized_(false),
  compressed_(false),
  stats_(),
  y_(),
  covar_()
  { }

  SumStatsLibrary(const PolymorphismData& dataset):
  numPops_(0),
  order_(0),
  initialized_(false),
  compressed_(false),
  stats_(),
  y_(),
  covar_()
  {
    init(dataset);
  }

  SumStatsLibrary(const OptionsContainer& opt):
  numPops_(opt.getNumberOfPopulations()),
  order_(opt.getOrder()),
  initialized_(false),
  compressed_(false),
  stats_(),
  y_(),
  covar_()
  {
    PolymorphismData dataset;
    dataset.parse(opt.getDataPath());
    dataset.computeSumStats();

    init(dataset);
  }

public:
  size_t getNumPops()
  {
    return numPops_;
  }

  size_t getNumPops() const
  {
    return numPops_;
  }

  size_t getOrder()
  {
    return order_;
  }

  size_t getOrder() const
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

  size_t getNumStats() const
  {
    return stats_.size();
  }

  const std::map<std::string, double>& getStats() const
  {
    return stats_;
  }

  double getStat(const std::string& name)
  {
    return stats_.at(name);
  }

  const Eigen::VectorXd& getYvec()
  {
    return y_;
  }

  const Eigen::MatrixXd getCovarMatrix()
  {
    return covar_;
  }

  void init();

  std::vector<std::string> splitString(const std::string& target, const std::string& query) const
  {
    std::vector<std::string> ret(0);
    boost::split(ret, target, boost::is_any_of(query));

    return ret;
  }

  size_t countInstances(const std::string& target, const std::string& query) const
  {
    std::string::difference_type count = std::count(std::begin(target), std::end(target), *query.c_str());
    return static_cast<size_t>(count);
  }

  size_t indexLookup(const std::string& moment) const
  {
    return std::distance(std::begin(stats_), stats_.find(moment));
  }

  void init(const PolymorphismData& dataset);

private:
  void includeHetStats_();

  void includeLdStats_();

  void compress_(); // exploits symmetry among statistics to reduce dimension

};

#endif

