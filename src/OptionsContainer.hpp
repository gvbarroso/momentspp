/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 14/10/2025
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
  std::string dataFilePath_;      // observed sum stats, for most recent Epoch
  std::string initStatsFilePath_; // e.g. steady-state sum stats for deep-most Epoch
  std::string numericalOptimizer_;
  std::string steadyState_; // eigen or pseudo, see Epoch class

  double toleranceOptim_; // for numerical optimization

  bool aliasOverEpochs_; // whether to alias parameters (r_*, u_*, s_*) over Epochs, see Model::compressParameters()
  bool aliasOverPops_;   // whether to alias parameters (r_*, u_*, s_*) over Populations, see Model::compressParameters()
  bool compressMoments_; // see moments_ vs compressedBasis_ inside SumStatsLibrary::initMoments_()
  bool computeCI_;
  bool verbose_;
  bool continuousTime_; // integration of expected sum. stats. (see Epoch and Model classes)

  size_t digits_; // if > 16 (DEFAULT), uses mpfr::mpreal, else uses double throughout the program execution
  size_t numThreads_;
  size_t timeSteps_; // to print intermediate values for Hl_*_* and Hr_*_* (see interval arg in Epoch::printMomentsIntermediate())
  std::vector<size_t> factorOrder_; // how many (1-2p) factors to include (one value per Epoch; if only one value is provided, it will be used for every Epoch)

  // for continuous-time integration, see Epoch class
  double dt_; // default time step for fixed integration
  double totalTimeIntegration_; // default total time for adaptive integration
  double toleranceIntegration_; // default error tolerance in adaptive (continuous-time) integration

  std::vector<std::string> momNames_; // moments to print at intermediate time steps

public:
  OptionsContainer(const std::map<std::string, std::string>& options):
  label_(bpp::ApplicationTools::getStringParameter("label", options, "moments++", "", 0, 4)),
  demesFilePath_(bpp::ApplicationTools::getAFilePath("demes_file", options, 0, 0, "", 0, "none", 0)),
  dataFilePath_(bpp::ApplicationTools::getAFilePath("obs_stats_file", options, false, true, "", false, "none", 4)),
  initStatsFilePath_(bpp::ApplicationTools::getAFilePath("init_stats_file", options, false, true, "", false, "none", 4)),
  numericalOptimizer_(bpp::ApplicationTools::getStringParameter("optimizer", options, "NewtonRhapson", "", true, 4)),
  steadyState_(bpp::ApplicationTools::getStringParameter("steady_state", options, "eigen", "", true, 4)),
  toleranceOptim_(bpp::ApplicationTools::getDoubleParameter("toleranceOptim", options, 1e-6, "", 0, 4)),
  aliasOverEpochs_(bpp::ApplicationTools::getParameter<bool>("alias_epochs_params", options, true, "", true, 4)),
  aliasOverPops_(bpp::ApplicationTools::getParameter<bool>("alias_pops_params", options, true, "", true, 4)),
  compressMoments_(bpp::ApplicationTools::getParameter<bool>("compress_moments", options, true, "", true, 4)),
  computeCI_(bpp::ApplicationTools::getParameter<bool>("ci", options, true, "", true, 4)),
  verbose_(bpp::ApplicationTools::getParameter<bool>("verbose", options, false, "", true, 4)),
  continuousTime_(bpp::ApplicationTools::getParameter<bool>("continuous_time", options, false, "", true, 4)),
  digits_(bpp::ApplicationTools::getParameter<size_t>("digits", options, 16, "", true, 4)),
  numThreads_(bpp::ApplicationTools::getParameter<size_t>("num_threads", options, std::thread::hardware_concurrency() / 2, "", true, 4)),
  timeSteps_(bpp::ApplicationTools::getParameter<size_t>("time_steps", options, 0, "", true, 4)),
  factorOrder_(bpp::ApplicationTools::getVectorParameter<size_t>("factor_order", options, ',', "10", "", true, 0)),
  dt_(bpp::ApplicationTools::getDoubleParameter("dt", options, 1e-3, "", true, 4)),
  totalTimeIntegration_(bpp::ApplicationTools::getDoubleParameter("time_integration", options, 1., "", 0, 4)),
  toleranceIntegration_(bpp::ApplicationTools::getDoubleParameter("tolerance_integration", options, 1e-6, "", 0, 4)),
  momNames_(bpp::ApplicationTools::getVectorParameter<std::string>("moms_intermediate", options, ',', "Hl_0_0,Hr_0_0", "", true, 4)) // NOTE check passing default values
  {
    if(label_ == "moments++")
      label_ = demesFilePath_.substr(0, demesFilePath_.find(".yaml")); // convenience

    if(digits_ < 16)
      std::cout << "WARNING: specified number of precision digits < 16. Will use double precision instead.\n\n";
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

  const std::string& getSteadyStateMethod() const
  {
    return steadyState_;
  }

  double getToleranceOptim() const
  {
    return toleranceOptim_;
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

  bool continuousTime() const
  {
    return continuousTime_;
  }

  bool highPrecision() const
  {
    return digits_ > 16; // default = 16 (double)
  }

  size_t getDigits() const
  {
    return digits_;
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

  double getDt() const
  {
    return dt_;
  }

  double getTotalTimeIntegration()
  {
    return totalTimeIntegration_;
  }

  double getToleranceIntegration()
  {
    return toleranceIntegration_;
  }

  const std::vector<std::string>& getMomNamesIntermediate() const
  {
    return momNames_;
  }
};

#endif
