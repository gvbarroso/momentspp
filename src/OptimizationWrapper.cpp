/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 29/07/2022
 *
 */


#include <iostream>

#include "SmcOptimizationWrapper.h"
#include "BackupListenerOv.h"
#include "Global.h"

using namespace std;
using namespace bpp;


std::shared_ptr<Model> OptimizationWrapper::selectBestModel()
{
  size_t bestModelIndex = 0;   
  double bestAic = listOfModels_[0]->getAic();

  for(size_t i = 1; i < listOfModels_.size(); ++i)
  {
    if(listOfModels_[i]->getAic() < bestAic)
    {
      bestAic = listOfModels_[i]->getAic();
      bestModelIndex = i;
    }
  }

  return listOfModels_[bestModelIndex];
}
  

void OptimizationWrapper::optimizeParameters() {

  bpp::ParameterList params(bestParameters_);
  
  createAndFitModels_(params);
  //updates params
  std::shared_ptr<SplinesModel> bestModel = selectBestModel();

  params.matchParametersValues(mmsmc_->getParameters());
  params.matchParametersValues(mmsmc_->getLambdaVector());

  for(size_t i = 0; i < mmsmc_->getParameterScalings().size(); ++i) {
    params.matchParametersValues(mmsmc_->getParameterScalings()[i]->getIndependentParameters());
  }
  
  for(size_t i = 0; i < mmsmc_->getParameterTransitions().size(); ++i) {
    params.matchParametersValues(mmsmc_->getParameterTransitions()[i]->getParameters());
  }
  
  fireUpdateBestValues_(bestSplines.get(), params);
}

void OptimizationWrapper::writeEstimatesToFile(const bpp::ParameterList& params, double aic) {
    
  //writes params in a simple txt format that is easy to parse
  std::ofstream parameterEstimates;
  
  parameterEstimates.open("estimates.txt");
  parameterEstimates << "AIC = " << aic << endl << endl;
  
  std::vector<std::string> parameterNames = params.getParameterNames();
  
  for(size_t i = 0; i < params.size(); ++i)
  {
    parameterEstimates << params.getParameterNames()[i];
    parameterEstimates << " ";

    parameterEstimates << params.getParameterValue(params.getParameterNames()[i]);
    parameterEstimates << std::endl;
  }
  
  parameterEstimates.close();    
}

void OptimizationWrapper::fireUpdateBestValues_(Model* bestModel, const bpp::ParameterList& params) {

  if(bestModel->getAic() < bestAic_)
  {
    bestAic_ = bestSm->getAic();
    bestParameters_.matchParametersValues(params);
    
    bpp::ParameterList estimates = bestModel->getParameters();
    writeEstimatesToFile(estimates, bestModel->getAic());
  }
  
  else
    std::cout << "WARNING!!! Optimisation did not reduce AIC!\n";
}  

void OptimizationWrapper::createAndFitModels_(ParameterList& params) {

  unsigned int minNumEpochs = options_->getMinNumberOfEpochs();
  unsigned int maxNumEpochs = options_->getMaxNumberOfEpochs();

  for(size_t i = minNumEpochs; i <= maxNumEpochs; ++i)
  {
    // adds N(t) parameters based on i
    addNeParams(params); // TODO
    std::shared_ptr<Model> model(params);

    fitModel_(model.get());
    listOfModels_.push_back(model);
  }
}
    
void OptimizationWrapper::fitModel_(Model* model) {
  
  std::cout << "\nOptimizing model with the following parameters:\n";

  model->fetchModelParameters().printParameters(std::cout);

  bpp::ReparametrizationFunctionWrapper rfw(model, model->fetchModelParameters());
  bpp::ThreePointsNumericalDerivative tpnd(model);
   
  std::unique_ptr<bpp::Optimizer> chosenOptimizer;

  if(options_->getOptimizer() == "Powell")
    chosenOptimizer.reset(new PowellMultiDimensions(&rfw));

  else if(options_->getOptimizer() == "NewtonRhapson")
  {
    tpnd.setParametersToDerivate(model->fetchNonSplinesParameters().getParameterNames());
    
    tpnd.enableFirstOrderDerivatives(true);
    tpnd.enableSecondOrderDerivatives(true);
    tpnd.enableSecondOrderCrossDerivatives(false); // Pseudo-Newton
    
    chosenOptimizer.reset(new bpp::PseudoNewtonOptimizer(&tpnd));
  }

  std::unique_ptr<StlOutputStream> profiler;
  std::unique_ptr<StlOutputStream> messenger;

  std::string optimProfile = "optim_profile.txt";
  profiler.reset(new StlOutputStream(new std::ofstream(optimProfile, ios::out)));
  chosenOptimizer->setProfiler(profiler.get());
  
  std::string optimMsgs = "optim_messages.txt";
  messenger.reset(new StlOutputStream(new std::ofstream(optimMsgs, ios::out)));
  chosenOptimizer->setMessageHandler(messenger.get());
    
  if(options_->getOptimizer() == "Powell")
    chosenOptimizer->init(rfw.getParameters());
  
  else if(options_->getOptimizer() == "NewtonRhapson")
  {
    chosenOptimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
    chosenOptimizer->init(model->fetchModelParameters());
  }
  
  std::unique_ptr<bpp::FunctionStopCondition> stopCond;
  stopCond.reset(new bpp::FunctionStopCondition(chosenOptimizer.get(), options_->getFunctionTolerance()));
  
  chosenOptimizer->setStopCondition(*stopCond);
  
  // handling optimisation issues, e.g. Powell: line minimization failing
  try
  {
    BackupListenerOv blo("backup_params.txt");
    chosenOptimizer->addOptimizationListener(&blo);
    chosenOptimizer->optimize();
  }

  catch(bpp::Exception& e)
  {
    std::cout << "\nWarning!! Error during optimization, convergence might not be reached.\n";
    std::cout << "moments++ will proceed, but check log files.\n";
  }
  
  model->computeAic();

  std::cout << "\n\nAIC = " << std::setprecision(6) << model->getAic() << std::endl;

  bpp::ParameterList optimParams(model->getParameters());
  std::cout << "\nOptimized parameters:\n";
  optimParams.printParameters(std::cout);
}


