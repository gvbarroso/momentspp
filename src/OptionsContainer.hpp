/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 04/04/2023
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


class OptionsContainer
{

private:
  std::string label_;
  std::string demesFilePath_;
  std::string dataFilePath_; // or observed sum stats
  std::string numericalOptimizer_;

  std::vector<double> initMij_;
  std::vector<double> initPopSizes_;

  double initMu_;
  double initR_;
  double tolerance_; // for numerical optimization

  bool compressMoments_; // see moments_ vs compressedBasis_ inside SumStatsLibrary::initMoments_()
  bool computeCI_;

  size_t order_; // sample order of summary statistics
  size_t numEpochs_;
  size_t numPops_;
  size_t totalNumberOfGenerations_;
  size_t numThreads_;

public:
  OptionsContainer(const std::map<std::string, std::string>& options):
  label_(bpp::ApplicationTools::getStringParameter("label", options, "moments++", "", 0, 4)),
  demesFilePath_(bpp::ApplicationTools::getAFilePath("demes_file", options, 0, 0, "", 0, "none", 0)),
  dataFilePath_(bpp::ApplicationTools::getAFilePath("stats_file", options, false, true, "", false, "none", 0)),
  numericalOptimizer_(bpp::ApplicationTools::getStringParameter("optimizer", options, "Powell", "", true, 4)),
  initMij_(bpp::ApplicationTools::getVectorParameter<double>("mij", options, ',', "none")),
  initPopSizes_(bpp::ApplicationTools::getVectorParameter<double>("Ni", options, ',', "none")),
  initMu_(bpp::ApplicationTools::getDoubleParameter("mu", options, 1e-8)),
  initR_(bpp::ApplicationTools::getDoubleParameter("r", options, 1e-8)),
  tolerance_(bpp::ApplicationTools::getDoubleParameter("tolerance", options, 1e-6, "", 0, 4)),
  compressMoments_(bpp::ApplicationTools::getParameter<bool>("compress_moments", options, true, "", true, 0)),
  computeCI_(bpp::ApplicationTools::getParameter<bool>("ci", options, true, "", true, 4)),
  order_(bpp::ApplicationTools::getParameter<size_t>("order", options, 2, "", true, 4)),
  numEpochs_(bpp::ApplicationTools::getParameter<size_t>("num_epochs", options, 1)),
  numPops_(bpp::ApplicationTools::getParameter<size_t>("num_pops", options, 1)),
  totalNumberOfGenerations_(bpp::ApplicationTools::getParameter<size_t>("total_gen", options, 1)),
  numThreads_(bpp::ApplicationTools::getParameter<size_t>("num_threads", options, std::thread::hardware_concurrency() / 2, "", true, 4))
  {
    if(numPops_ != initPopSizes_.size())
      throw bpp::Exception("OptionsContainer::num_pops does not match length of Ni parameters!");

    if(numPops_ > 1)
      if(numPops_ * (numPops_ - 1) != initMij_.size())
        throw bpp::Exception("OptionsContainer::num_pops is not compatible with length of mij parameters!");
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

  const std::string& getOptimMethod() const
  {
    return numericalOptimizer_;
  }

  const std::vector<double>& getInitMig() const
  {
    return initMij_;
  }

  const std::vector<double>& getInitPopSizes() const
  {
    return initPopSizes_;
  }

  double getInitMu() const
  {
    return initMu_;
  }

  double getInitR() const
  {
    return initR_;
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

  size_t getOrder() const
  {
    return order_;
  }

  size_t getNumEpochs() const
  {
    return numEpochs_;
  }

  size_t getNumPops() const
  {
    return numPops_;
  }

  size_t getTotalNumberOfGenerations() const
  {
    return totalNumberOfGenerations_;
  }

  size_t getNumThreads() const
  {
    return numThreads_;
  }

};

#endif
