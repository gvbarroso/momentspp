/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 08/06/2023
 *
 */


#ifndef _OPTIONSCONTAINER_H_
#define _OPTIONSCONTAINER_H_

#include <string>
#include <vector>
#include <map>
#include <limits>
#include <thread>

#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Text/TextTools.h>

class OptionsContainer
{

private:
  std::string label_;
  std::string demesFilePath_;
  std::string dataFilePath_; // or observed sum stats
  std::string numericalOptimizer_;

  double tolerance_; // for numerical optimization
  bool compressMoments_; // see moments_ vs compressedBasis_ inside SumStatsLibrary::initMoments_()
  bool computeCI_;
  size_t numThreads_;
  size_t factorOrder_; // how many (1-2p) factors to include

public:
  OptionsContainer(const std::map<std::string, std::string>& options):
  label_(bpp::ApplicationTools::getStringParameter("label", options, "moments++", "", 0, 4)),
  demesFilePath_(bpp::ApplicationTools::getAFilePath("demes_file", options, 0, 0, "", 0, "none", 0)),
  dataFilePath_(bpp::ApplicationTools::getAFilePath("stats_file", options, false, true, "", false, "none", 0)),
  numericalOptimizer_(bpp::ApplicationTools::getStringParameter("optimizer", options, "NewtonRhapson", "", true, 4)),
  tolerance_(bpp::ApplicationTools::getDoubleParameter("tolerance", options, 1e-6, "", 0, 4)),
  compressMoments_(bpp::ApplicationTools::getParameter<bool>("compress_moments", options, true, "", true, 0)),
  computeCI_(bpp::ApplicationTools::getParameter<bool>("ci", options, true, "", true, 4)),
  numThreads_(bpp::ApplicationTools::getParameter<size_t>("num_threads", options, std::thread::hardware_concurrency() / 2, "", true, 4)),
  factorOrder_(bpp::ApplicationTools::getParameter<size_t>("factor_order", options, 1, "", true, 0))
  { }
  
public:
  const std::string& getLabel() const
  {
    return label_;
  }

  const std::string& getDemesFilePath() const
  {
    return demesFilePath_;
  }

  const std::string& getDataFilePath() const
  {
    return dataFilePath_;
  }

  const std::string& getOptimMethod() const
  {
    return numericalOptimizer_;
  }

  double getTolerance() const
  {
    return tolerance_;
  }

  bool compressMoments() const
  {
    return compressMoments_;
  }

  bool computeCI() const
  {
    return computeCI_;
  }

  size_t getNumThreads() const
  {
    return numThreads_;
  }

  size_t getFactorOrder() const
  {
    return factorOrder_;
  }

};

#endif
