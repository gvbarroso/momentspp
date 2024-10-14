/*
 * Authors: Gustavo V. Barroso
 * Created: 06/09/2022
 * Last modified: 30/05/2024
 *
 */

#ifndef _DATA_H_
#define _DATA_H_

#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <memory>

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>

#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/zlib.hpp>


#include "Population.hpp"
#include "SumStatsLibrary.hpp"
#include "OptionsContainer.hpp"

class Data // observed data
{

private:
  SumStatsLibrary ssl_; // the moments in Data naturally refer to Population indices among sampled individuals
  Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> covar_; // covariance matrix of observed sum stats (from sampled populations); ssl_.fetchYvec() will get the expectations
  std::map<std::string, double> variances_; // bootstrapped, moment name -> moment var

public:
  Data(const std::string& file):
  ssl_(),
  covar_(),
  variances_()
  {
    parse_(file);
  }

public:
  const SumStatsLibrary& getSumStatsLibrary()
  {
    return ssl_;
  }

  const Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic>& getCovarMatrix()
  {
    return covar_;
  }

  Eigen::Matrix<long double, Eigen::Dynamic, 1> getY()
  {
    return ssl_.fetchYvec();
  }

private:
  void parse_(const std::string& file);

};

#endif
