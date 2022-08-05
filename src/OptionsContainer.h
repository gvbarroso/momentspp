/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 29/07/2022
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

  double tolerance_;

  size_t numberOfPopulations_;
  size_t minNumberOfEpochs_;
  size_t maxNumberOfEpochs_;
  size_t numberOfThreads_;
  
public:
  OptionsContainer(std::map<std::string, std::string> options):
  dataPath_(bpp::ApplicationTools::getAFilePath("data_path", options, "none")),
  numericalOptimizer_(bpp::ApplicationTools::getStringParameter("optimizer", options, "Powell", "", true, 4)),
  tolerance_(bpp::ApplicationTools::getDoubleParameter("tolerance", options, 1e-6)),
  numberOfPopulations_(bpp::ApplicationTools::getParameter<size_t>("num_pops", options, "none")),
  minNumberOfEpochs_(bpp::ApplicationTools::getParameter<size_t>("min_epochs", options, "none")),
  maxNumberOfEpochs_(bpp::ApplicationTools::getParameter<size_t>("max_epochs", options, "none")),
  numberOfThreads_(bpp::ApplicationTools::getParameter<size_t>("number_threads", options,
                                                               std::thread::hardware_concurrency(),
                                                               "", true, 4)),
  { }
  
public:
  const std::string& getDataPath() const
  {
    return dataPath_;
  }

  const std::string& getOptimizer() const
  {
    return numericalOptimizer_;
  }

  size_t getMinNumberOfEpochs()
  {
    return initNumberOfKnots_;
  }
  
  size_t getMaxNumberOfEpochs()
  {
    return maxNumberOfKnots_;
  }

  size_t getNumberOfThreads()
  {
    return numberOfThreads_;
  }

  double getTolerance()
  {
    return tolerance_;
  }

};

#endif
