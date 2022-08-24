/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 24/08/2022
 *
 */


#ifndef _OPTIMIZATIONWRAPPER_H_
#define _OPTIMIZATIONWRAPPER_H_

#include <string>
#include <vector>
#include <limits>

#include <Bpp/Numeric/Function/ReparametrizationFunctionWrapper.h>
#include <Bpp/Numeric/Function/PowellMultiDimensions.h>
#include <Bpp/Numeric/Function/Bfgs.h>
#include <Bpp/Numeric/Function/ThreePointsNumericalDerivative.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Numeric/ParameterList.h>
#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Phyl/PseudoNewtonOptimizer.h>

#include "Model.hpp"
#include "OptionsContainer.hpp"


class OptimizationWrapper
{
    
private:
  OptionsContainer options_;
  std::vector<std::shared_ptr<Model>> listOfModels_;

  bpp::ParameterList bestParameters_;
  double bestAic_;
  
public:
  OptimizationWrapper(const OptionsContainer& options):
  options_(options),
  listOfModels_(0),
  bestParameters_(),
  bestAic_(std::numeric_limits<double>::max()) 
  {

  }
  
public:
    
  const bpp::ParameterList& getBestParameters() const
  {
    return bestParameters_;
  }
    
  std::vector<std::shared_ptr<Model>>& getListOfModels()
  {
    return listOfModels_;
  }

  const std::vector<std::shared_ptr<Model>>& getListOfModels() const
  {
    return listOfModels_;
  }
  
  double getAic()
  {
    return bestAic_;
  }
  
  std::shared_ptr<Model> selectBestModel();
  
  void optimize();
  
  //to resume optimisation after a crash:
  void optimizeParameters(const bpp::ParameterList& backupParams);
  
  void writeEstimatesToFile(std::shared_ptr<Model> model);

private:
  void fireUpdateBestValues_(Model* bestModel, const bpp::ParameterList& params);
  
  void createAndFitModels_(bpp::ParameterList& nonSplinesParams);
    
  void fitModel_(Model* model);
  
};

#endif
