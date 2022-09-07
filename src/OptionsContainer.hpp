/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 06/09/2022
 *
 */


#ifndef _OPTIONSCONTAINER_H_
#define _OPTIONSCONTAINER_H_

#include <string>
#include <vector>
#include <limits>
#include <thread>

#include <Bpp/App/ApplicationTools.h>


class OptionsContainer
{

private:
  std::string dataPath_;
  std::string numericalOptimizer_;

  double tolerance_; // for numerical optimization

  bool computeCI_;
  bool resume_;

  size_t order_; // of summary statistics
  size_t numberOfPopulations_;
  size_t minNumberOfEpochs_;
  size_t maxNumberOfEpochs_;
  size_t numberOfThreads_;
  
public:
  OptionsContainer(const std::map<std::string, std::string>& options):
  dataPath_(bpp::ApplicationTools::getAFilePath("data_path", options, "none")),
  numericalOptimizer_(bpp::ApplicationTools::getStringParameter("optimizer", options, "Powell", "", true, 4)),
  tolerance_(bpp::ApplicationTools::getDoubleParameter("tolerance", options, 1e-6)),
  computeCI_(bpp::ApplicationTools::getParameter<bool>("ci", options, true)),
  resume_(bpp::ApplicationTools::getParameter<bool>("resume", options, false)),
  order_(bpp::ApplicationTools::getParameter<size_t>("order", options, 2)),
  numberOfPopulations_(bpp::ApplicationTools::getParameter<size_t>("num_pops", options, 1)),
  minNumberOfEpochs_(bpp::ApplicationTools::getParameter<size_t>("min_epochs", options, 1)),
  maxNumberOfEpochs_(bpp::ApplicationTools::getParameter<size_t>("max_epochs", options, 1)),
  numberOfThreads_(bpp::ApplicationTools::getParameter<size_t>("number_threads", options,
                                                               std::thread::hardware_concurrency(),
                                                               "", true, 4))
  { }
  
public:
  const std::string& getDataPath() const
  {
    return dataPath_;
  }

  const std::string& getOptimMethod() const
  {
    return numericalOptimizer_;
  }

  double getTolerance() const
  {
    return tolerance_;
  }

  bool computeCI() const
  {
    return computeCI_;
  }

  bool resume() const
  {
    return resume_;
  }

  size_t getOrder() const
  {
    return order_;
  }

  size_t getNumberOfPopulations() const
  {
    return numberOfPopulations_;
  }

  size_t getMinNumberOfEpochs() const
  {
    return minNumberOfEpochs_;
  }
  
  size_t getMaxNumberOfEpochs() const
  {
    return maxNumberOfEpochs_;
  }

  size_t getNumberOfThreads() const
  {
    return numberOfThreads_;
  }

};

#endif
