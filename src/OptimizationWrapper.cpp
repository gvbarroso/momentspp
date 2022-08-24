/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 24/08/2022
 *
 */


#include <iostream>

#include "OptimizationWrapper.h"
#include "BackupListenerOv.h"

std::shared_ptr<Model> OptimizationWrapper::selectBestModel()
{
  size_t index = 0;
  double bestAic = listOfModels_[0]->getAic();

  for(size_t i = 1; i < listOfModels_.size(); ++i)
  {
    if(listOfModels_[i]->getAic() < bestAic)
    {
      bestAic = listOfModels_[i]->getAic();
      index = i;
    }
  }

  return listOfModels_[index];
}

void OptimizationWrapper::optimize()
{

  bpp::ParameterList params(bestParameters_);
  
  createAndFitModels_(params);
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

void OptimizationWrapper::writeEstimatesToFile(std::shared_ptr<Model> model) {
    
  //writes params in a simple txt format that is easy to parse
  std::ofstream file;
  
  file.open("moments++estimates.txt");
  file << "AIC = " << model->aic() << std::endl << std::endl;

  for(size_t i = 0; i < model->getParameters().size(); ++i)
    file << getParameters()[i].getName() << "\t" << getParameters()[i].getValue() << "\n";
  
  file.close();
}

void OptimizationWrapper::fireUpdateBestValues_(Model* model) {

  if(model->aic() < bestAic_)
  {
    bestAic_ = model->aic();
    bestParameters_ = model->getParameters();
    
    writeEstimatesToFile(model);
  }
}  

void OptimizationWrapper::createAndFitModels_(ParameterList& params) {

  unsigned int minNumEpochs = options_.getMinNumberOfEpochs();
  unsigned int maxNumEpochs = options_.getMaxNumberOfEpochs();

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

  model->getIndependentParameters().printParameters(std::cout);

  bpp::Optimizer chosenOptimizer;

  std::unique_ptr<bpp::StlOutputStream> profiler;
  std::unique_ptr<bpp::StlOutputStream> messenger;

  std::string optimProfile = "profile.txt";
  profiler.reset(new bpp::StlOutputStream(new std::ofstream(optimProfile, std::ios::out)));
  chosenOptimizer.setProfiler(profiler.get());
  
  std::string optimMsgs = "messages.txt";
  messenger.reset(new StlOutputStream(new std::ofstream(optimMsgs, std::ios::out)));
  chosenOptimizer.setMessageHandler(messenger.get());
    
  if(options_.getOptimizer() == "Powell")
  {
    bpp::ReparametrizationFunctionWrapper rfw(model, model->getIndependentParameters());

    chosenOptimizer.reset(new PowellMultiDimensions(&rfw));
    chosenOptimizer.init(rfw.getParameters());
  }

  else if(options_.getOptimizer() == "NewtonRhapson")
  {
    bpp::ThreePointsNumericalDerivative tpnd(model);

    tpnd.setParametersToDerivate(model->fetchNonSplinesParameters().getParameterNames());
    tpnd.enableFirstOrderDerivatives(true);
    tpnd.enableSecondOrderDerivatives(true);
    tpnd.enableSecondOrderCrossDerivatives(false); // Pseudo-Newton

    chosenOptimizer.reset(new bpp::PseudoNewtonOptimizer(&tpnd));
    chosenOptimizer.setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
    chosenOptimizer.init(model->fetchModelParameters());
  }
  
  bpp::FunctionStopCondition stopCond;
  stopCond.reset(new bpp::FunctionStopCondition(chosenOptimizer.get(), options_.getFunctionTolerance()));
  
  chosenOptimizer.setStopCondition(*stopCond);
  
  // handling optimisation issues
  try
  {
    BackupListenerOv blo("backup_params.txt");
    chosenOptimizer.addOptimizationListener(&blo);
    chosenOptimizer.optimize();
  }

  catch(bpp::Exception& e)
  {
    std::cout << "\nWarning!! Error during optimization, convergence might not be reached.\n";
    std::cout << "moments++ will proceed, but check log files.\n";
  }
  
  model->computeAic();

  std::cout << "\n\nAIC = " << std::setprecision(6) << model->aic() << std::endl;

  std::cout << "\nOptimized parameters:\n";
  model->getParameters().printParameters(std::cout);
}


