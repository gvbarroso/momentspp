/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 28/09/2022
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
  std::string popsFilePath_;
  std::string dataFilePath_;
  std::string numericalOptimizer_;

  double tolerance_; // for numerical optimization

  bool computeCI_;

  size_t order_; // of summary statistics
  size_t totalNumberOfGenerations_;
  size_t numberOfThreads_;

public:
  OptionsContainer(const std::map<std::string, std::string>& options):
  popsFilePath_(bpp::ApplicationTools::getAFilePath("pop_file", options, "none")),
  dataFilePath_(bpp::ApplicationTools::getAFilePath("data_file", options, "none")),
  numericalOptimizer_(bpp::ApplicationTools::getStringParameter("optimizer", options, "Powell", "", true, 4)),
  tolerance_(bpp::ApplicationTools::getDoubleParameter("tolerance", options, 1e-6)),
  computeCI_(bpp::ApplicationTools::getParameter<bool>("ci", options, true)),
  order_(bpp::ApplicationTools::getParameter<size_t>("order", options, 2)),
  totalNumberOfGenerations_(bpp::ApplicationTools::getParameter<size_t>("total_gen", options, 1)),
  numberOfThreads_(bpp::ApplicationTools::getParameter<size_t>("number_threads", options,
                                                               std::thread::hardware_concurrency(),
                                                               "", true, 4))
  { }
  
public:
  const std::string& getPopsFilePath() const
  {
    return popsFilePath_;
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

  bool computeCI() const
  {
    return computeCI_;
  }

  size_t getOrder() const
  {
    return order_;
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
