/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 25/08/2022
 *
 */


#include <iostream>

#include "SumStatsLibrary.hpp"
#include "OptimizationWrapper.hpp"
#include "BackupListenerOv.hpp"


void OptimizationWrapper::optimize(const SumStatsLibrary& sslib)
{
  std::string name = "";

  unsigned int minNumEpochs = options_.getMinNumberOfEpochs();
  unsigned int maxNumEpochs = options_.getMaxNumberOfEpochs();

  for(size_t i = minNumEpochs; i <= maxNumEpochs; ++i)
  {
    name = "moments++" + "_model_" + bpp::ToString(i);

    // build parameter list with constraints
    bpp::ParameterList params;

    // build operators NOTE mind correct order for multiplying matrices into combined operator
    std::vector<Operator*> operators(0);

    // e.g.
    Drift* driftOperator = new Drift(params, sslib);
    operators.push_back(driftOperator);

    Model* model = new Model(name, operators, params, sslib);

    /*
     * decides whether to freeze parameters
     */

    fitModel_(model.get());
    writeEstimatesToFile_(model);

    // clean-up
    delete model;

    for(auto it = std::begin(operators); it != std::end(operators); ++it)
    {
      delete *it;
      it = operators.erase(it);
    }
  }
}

void OptimizationWrapper::fitModel_(Model* model) {
  
  std::cout << "\nOptimizing model with the following parameters:\n";

  model->getUnfrozenParameters().printParameters(std::cout);

  bpp::Optimizer chosenOptimizer;

  bpp::StlOutputStream profiler;
  bpp::StlOutputStream messenger;

  std::string optimProfile = "profile.txt";
  profiler.reset(new bpp::StlOutputStream(new std::ofstream(optimProfile, std::ios::out)));
  chosenOptimizer.setProfiler(&profiler);
  
  std::string optimMsgs = "messages.txt";
  messenger.reset(new StlOutputStream(new std::ofstream(optimMsgs, std::ios::out)));
  chosenOptimizer.setMessageHandler(&messenger);
    
  if(options_.getOptimMethod() == "Powell")
  {
    bpp::ReparametrizationFunctionWrapper rfw(model, model->getUnfrozenParameters());

    chosenOptimizer.reset(new PowellMultiDimensions(&rfw));
    chosenOptimizer.init(rfw.getParameters());
  }

  else if(options_.getOptimMethod() == "NewtonRhapson")
  {
    bpp::ThreePointsNumericalDerivative tpnd(model);

    tpnd.setParametersToDerivate(model->getUnfrozenParameters().getParameterNames());
    tpnd.enableFirstOrderDerivatives(true);
    tpnd.enableSecondOrderDerivatives(true);
    tpnd.enableSecondOrderCrossDerivatives(false); // Pseudo-Newton

    chosenOptimizer.reset(new bpp::PseudoNewtonOptimizer(&tpnd));
    chosenOptimizer.setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
    chosenOptimizer.init(model->fetchModelParameters());
  }
  
  else
    throw bpp::Excepetion("OptimizationWrapper::Mis-specified numerical optimizer");

  bpp::FunctionStopCondition* stopCond =
  stopCond.reset(new bpp::FunctionStopCondition(&chosenOptimizer, options_.getFunctionTolerance()));
  
  chosenOptimizer.setStopCondition(*stopCond);
  
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
}

void OptimizationWrapper::writeEstimatesToFile_(Model* model)
{
  std::ofstream file;
  file.open(model->getName() + "_estimates.txt");

  file << "AIC = " << model->aic() << std::endl << std::endl;
  model->getParameterList().printParameters(file);

  file.close();
}
