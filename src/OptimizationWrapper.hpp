/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 29/07/2022
 *
 */


#ifndef _OPTIMIZATIONWRAPPER_H_
#define _OPTIMIZATIONWRAPPER_H_

#include <string>
#include <vector>
#include <limits>

#include <Bpp/Numeric/Function/ReparametrizationFunctionWrapper.h>
#include <Bpp/Numeric/Function/PowellMultiDimensions.h>
#include <Bpp/Numeric/Function/ThreePointsNumericalDerivative.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/ParameterList.h>
#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Phyl/PseudoNewtonOptimizer.h>

#include "Model.h"
#include "OptionsContainer.h"


class OptimizationWrapper
{
    
private:
  std::shared_ptr<OptionsContainer> options_;
  std::vector<std::shared_ptr<Model>> listOfModels_;

  bpp::ParameterList bestParameters_;
  double bestAic_;
  
public:
  OptimizationWrapper(std::shared_ptr<OptionsContainer> options):
  options_(options),
  listOfModels_(0),
  bestParameters_(),
  bestAic_(std::numeric_limits<double>::max()) 
  {
    //standard SMC parameters
    bestParameters_.addParameters(mmsmc -> getParameters());
    bestParameters_.addParameters(mmsmc -> getLambdaVector());
    
    //add Markov-modulation parameters
    for(size_t i = 0; i < mmsmc -> getParameterScalings().size(); ++i) {
        
      //if hotspot model, only bring hotspot intensity to optimization (ie discard PMF probs.)
      if(mmsmc -> getParameterTransitions()[i] -> getHeterogeneousRateModel() == "Hotspot") {
        bestParameters_.addParameter(mmsmc -> getParameterScalings()[i] -> getParameter("V2"));  
      }
      
      else if(mmsmc -> getParameterTransitions()[i] -> getHeterogeneousRateModel() == "Gamma") {
          
        if(mmsmc -> getParameterScalings()[i] -> getNumberOfCategories() > 1) {
          bestParameters_.addParameters(mmsmc -> getParameterScalings()[i] -> getIndependentParameters()); 
        }
      }
      
      else if(mmsmc -> getParameterTransitions()[i] -> getHeterogeneousRateModel() == "Gamma+Hotspot") {
          
        if(mmsmc -> getParameterScalings()[i] -> getNumberOfCategories() > 1) {
          bestParameters_.addParameters(mmsmc -> getParameterScalings()[i] -> getIndependentParameters()); 
        }
        
        else {
          bestParameters_.addParameter(mmsmc -> getParameterScalings()[i] -> getParameter("heat"));
        }
      }
    }
    
    for(size_t i = 0; i < mmsmc -> getParameterTransitions().size(); ++i) {
        
      if(mmsmc -> getParameterTransitions()[i] -> getHeterogeneousRateModel() == "Gamma") {
          
        if(mmsmc -> getParameterScalings()[i] -> getNumberOfCategories() > 1) {
          bestParameters_.addParameters(mmsmc -> getParameterTransitions()[i] -> getParameters());
        }
      }
      
      else {
        bestParameters_.addParameters(mmsmc -> getParameterTransitions()[i] -> getParameters());
      }
    }
  }
  
public:
    
  bpp::ParameterList& getBestParameters()
  {
    return bestParameters_;
  }
  const bpp::ParameterList& getBestParameters() const
  {
    return bestParameters_;
  }
    
  std::vector<std::shared_ptr <Model>>& getListOfModels()
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

  void stepwiseExpectationMaximization();
  
  void optimizeParameters();
  
  //to resume optimisation after a problem:
  void optimizeParameters(const bpp::ParameterList& backupParams);
  
  void writeEstimatesToFile(const bpp::ParameterList& params, double AIC);

  void writeDemographyToFile();

private:
  void fireUpdateBestValues_(Model* bestModel, const bpp::ParameterList& params);
  
  void createAndFitSplinesModels_(bpp::ParameterList& nonSplinesParams);
    
  void fitModel_(Model* model);
  
};

#endif
