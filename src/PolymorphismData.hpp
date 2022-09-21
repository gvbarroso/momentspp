/*
 * Authors: Gustavo V. Barroso
 * Created: 06/09/2022
 * Last modified: 20/09/2022
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

#include Mo
#include "OptionsContainer.hpp"

class PolymorphismData {
private:
  size_t numPops_; // number of sampled populations in present time, set in parse() method
  size_t order_;

  std::map<std::string, double> expectations_;
  std::map<std::string, double> variances_; // bootstrapped

  Eigen::VectorXd y_; // the Eigen representation of the observed vector of summary statistics (from sampled populations)
  Eigen::MatrixXd covar_; // covariance matrix of observed sum stats (from sampled populations)

public:
  PolymorphismData():
  numPops_(1),
  order_(2),
  expectations_(),
  variances_(),
  y_(),
  covar_()
  { }

  PolymorphismData(const OptionsContainer& opt):
  numPops_(1),
  order_(opt.getOrder()),
  expectations_(),
  variances_(),
  y_(),
  covar_()
  { }

public:
  size_t getNumPops()
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
