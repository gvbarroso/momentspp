/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 18/06/2024
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
  std::string dataFilePath_; // observed sum stats, for most recent Epoch
  std::string initStatsFilePath_; // e.g. steady-state sum stats for deep-most Epoch
  std::string numericalOptimizer_;

  double tolerance_; // for numerical optimization

  bool aliasOverEpochs_; // whether to alias parameters (r_*, u_*, s_*) over Epochs, see Model::compressParameters()
  bool aliasOverPops_; // whether to alias parameters (r_*, u_*, s_*) over Populations, see Model::compressParameters()
  bool compressMoments_; // see moments_ vs compressedBasis_ inside SumStatsLibrary::initMoments_()
  bool computeCI_;
  bool verbose_;

  size_t numThreads_;
  size_t timeSteps_; // after how number of generations to print intermediate values for Hl_*_* and Hr_*_* (see interval arg in Epoch::printHetMomentsIntermediate())
  std::vector<size_t> factorOrder_; // how many (1-2p) factors to include (one value per Epoch; if only one value is provided, the same will be used)

public:
  OptionsContainer(const std::map<std::string, std::string>& options):
  label_(bpp::ApplicationTools::getStringParameter("label", options, "moments++", "", 0, 4)),
  demesFilePath_(bpp::ApplicationTools::getAFilePath("demes_file", options, 0, 0, "", 0, "none", 0)),
  dataFilePath_(bpp::ApplicationTools::getAFilePath("obs_stats_file", options, false, true, "", false, "none", 4)),
  initStatsFilePath_(bpp::ApplicationTools::getAFilePath("init_stats_file", options, false, true, "", false, "none", 4)),
  numericalOptimizer_(bpp::ApplicationTools::getStringParameter("optimizer", options, "NewtonRhapson", "", true, 4)),
  tolerance_(bpp::ApplicationTools::getDoubleParameter("tolerance", options, 1e-6, "", 0, 4)),
  aliasOverEpochs_(bpp::ApplicationTools::getParameter<bool>("alias_epochs_params", options, true, "", true, 4)),
  aliasOverPops_(bpp::ApplicationTools::getParameter<bool>("alias_pops_params", options, true, "", true, 4)),
  compressMoments_(bpp::ApplicationTools::getParameter<bool>("compress_moments", options, true, "", true, 4)),
  computeCI_(bpp::ApplicationTools::getParameter<bool>("ci", options, true, "", true, 4)),
  verbose_(bpp::ApplicationTools::getParameter<bool>("verbose", options, false, "", true, 4)),
  numThreads_(bpp::ApplicationTools::getParameter<size_t>("num_threads", options, std::thread::hardware_concurrency() / 2, "", true, 4)),
  timeSteps_(bpp::ApplicationTools::getParameter<size_t>("time_steps", options, 1000, "", true, 4)),
  factorOrder_(bpp::ApplicationTools::getVectorParameter<size_t>("factor_order", options, ',', "0", "", true, 0))
  {
    if(label_ == "moments++")
      label_ = demesFilePath_.substr(0, demesFilePath_.find(".yaml")); // convenience
  }
  
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

  const std::string& getInitStatsFilePath() const
  {
    return initStatsFilePath_;
  }

  const std::string& getOptimMethod() const
  {
    return numericalOptimizer_;
  }

  double getTolerance() const
  {
    return tolerance_;
  }

  bool aliasEpochsParams() const
  {
    return aliasOverEpochs_;
  }

  bool aliasPopsParams() const
  {
    return aliasOverPops_;
  }

  bool compressMoments() const
  {
    return compressMoments_;
  }

  bool computeCI() const
  {
    return computeCI_;
  }

  bool verbose() const
  {
    return verbose_;
  }

  size_t getNumThreads() const
  {
    return numThreads_;
  }

  size_t getTimeSteps() const
  {
    return timeSteps_;
  }

  const std::vector<size_t>& getFactorOrder() const
  {
    return factorOrder_;
  }


};

#endif
