/*
 * Authors: Gustavo V. Barroso
 * Created: 06/09/2022
 * Last modified: 11/09/2022
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

#include "OptionsContainer.hpp"

class PolymorphismData {
private:
  std::vector<size_t> numPops_;
  size_t order_;

  std::map<std::string, double> expectations_;
  std::map<std::string, double> variances_; // bootstrapped

public:
  PolymorphismData():
  numPops_(1),
  order_(2),
  expectations_(),
  variances_()
  { }

  PolymorphismData(const OptionsContainer& opt):
  numPops_(opt.getNumbersOfPopulations()),
  order_(opt.getOrder()),
  expectations_(),
  variances_()
  { }

public:
  std::vector<size_t> getNumPops()
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

  void parse(const std::string& file);

  void computeSumStats();

};

#endif
