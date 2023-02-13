/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 14/12/2022
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
  //std::string initParamsFilePath_; // for creating multiple models
  std::string numericalOptimizer_;

  std::vector<double> initMij_;
  std::vector<double> initPopSizes_;

  double initMu_;
  double initR_;
  double tolerance_; // for numerical optimization

  bool computeCI_;

  size_t order_; // of summary statistics
  size_t numEpochs_;
  size_t numPops_;
  size_t totalNumberOfGenerations_;
  size_t numberOfThreads_;

public:
  OptionsContainer(const std::map<std::string, std::string>& options):
  label_(bpp::ApplicationTools::getStringParameter("label", options, "my_model")),
  demesFilePath_(bpp::ApplicationTools::getAFilePath("demes_file", options, "none")),
  dataFilePath_(bpp::ApplicationTools::getAFilePath("stats_file", options, "none")),
  //initParamsFilePath_(bpp::ApplicationTools::getAFilePath("params_file", options, "none")),
  numericalOptimizer_(bpp::ApplicationTools::getStringParameter("optimizer", options, "Powell", "", true, 4)),
  initMij_(bpp::ApplicationTools::getVectorParameter<double>("mij", options, ',', "none")),
  initPopSizes_(bpp::ApplicationTools::getVectorParameter<double>("Ni", options, ',', "none")),
  initMu_(bpp::ApplicationTools::getDoubleParameter("mu", options, 1e-8)),
  initR_(bpp::ApplicationTools::getDoubleParameter("r", options, 1e-9)),
  tolerance_(bpp::ApplicationTools::getDoubleParameter("tolerance", options, 1e-6)),
  computeCI_(bpp::ApplicationTools::getParameter<bool>("ci", options, true)),
  order_(bpp::ApplicationTools::getParameter<size_t>("order", options, 2)),
  numEpochs_(bpp::ApplicationTools::getParameter<size_t>("num_epochs", options, 1)),
  numPops_(bpp::ApplicationTools::getParameter<size_t>("num_pops", options, 1)),
  totalNumberOfGenerations_(bpp::ApplicationTools::getParameter<size_t>("total_gen", options, 1)),
  numberOfThreads_(bpp::ApplicationTools::getParameter<size_t>("number_threads", options,
                                                               std::thread::hardware_concurrency(),
                                                               "", true, 4))
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
  /*
  const std::string& getParamsFilePath() const
  {
    return initParamsFilePath_;
  }*/

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

  size_t getNumberOfThreads() const
  {
    return numberOfThreads_;
  }

};

#endif
