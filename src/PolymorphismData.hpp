/*
 * Authors: Gustavo V. Barroso
 * Created: 06/09/2022
 * Last modified: 28/09/2022
 *
 */

#ifndef _POLYMORPHISMDATA_H_
#define _POLYMORPHISMDATA_H_

#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <memory>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/zlib.hpp>

/*#include <Bpp/Seq/Io/Fasta.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/SiteTools.h>*/

#include "Population.hpp"
#include "Moment.hpp"
#include "OptionsContainer.hpp"

class PolymorphismData {
private:
  size_t order_;

  std::vector<std::map<size_t, std::shared_ptr<Population>>> popMaps_; // pop-id -> pop, one per epoch

  std::map<std::string, double> expectations_;
  std::map<std::string, double> variances_; // bootstrapped

  Eigen::VectorXd y_; // the Eigen representation of the observed vector of summary statistics (from sampled populations)
  Eigen::MatrixXd covar_; // covariance matrix of observed sum stats (from sampled populations)

public:
  PolymorphismData():
  order_(2),
  popMaps_(),
  expectations_(),
  variances_(),
  y_(),
  covar_()
  { }

  PolymorphismData(const OptionsContainer& opt, const std::vector<std::map<size_t, std::shared_ptr<Population>>>& popMaps):
  order_(opt.getOrder()),
  popMaps_(popMaps),
  expectations_(),
  variances_(),
  y_(),
  covar_()
  { }

public:
  size_t getNumPops(size_t epochIndex)
  {
    return popMaps_[epochIndex].size();
  }

  size_t getNumSampledPops()
  {
    return popMaps_.back().size(); // TODO accomodate ancient samples
  }

  size_t getOrder()
  {
    return order_;
  }

  size_t getOrder() const
  {
    return order_;
  }

  const std::vector<std::map<size_t, std::shared_ptr<Population>>> getPopMaps() const
  {
    return popMaps_;
  }

  const Eigen::VectorXd& getYvec()
  {
    return y_;
  }

  const Eigen::MatrixXd& getCovarMatrix()
  {
    return covar_;
  }

  void parse(const std::string& file);

  void computeSumStats(); // from sampled populations

};

#endif
