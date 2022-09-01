/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 31/08/2022
 *
 */


#include <iostream>
#include <memory>

#include "SumStatsLibrary.hpp"
#include "OptimizationWrapper.hpp"
#include "BackupListenerOv.hpp"


void OptimizationWrapper::optimize(const SumStatsLibrary& sslib)
{
  size_t minNumEpochs = options_.getMinNumberOfEpochs();
  size_t maxNumEpochs = options_.getMaxNumberOfEpochs();

  for(size_t i = minNumEpochs; i <= maxNumEpochs; ++i)
  {
    std::string name = "moments++" + "_model_" + bpp::ToString(i);

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

    fitModel_(model);
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

void OptimizationWrapper::fitModel_(Model* model)
{
  std::cout << "\nOptimizing model with the following parameters:\n";
  model->getUnfrozenParameters().printParameters(std::cout);

  std::unique_ptr<bpp::Optimizer> optimizer;

  if(options_.getOptimMethod() == "Powell")
  {
    bpp::ReparametrizationFunctionWrapper rfw(model, model->getUnfrozenParameters());

    optimizer.reset(new PowellMultiDimensions(&rfw));
    optimizer->init(rfw.getParameters());
  }

  else if(options_.getOptimMethod() == "NewtonRhapson")
  {
    bpp::ThreePointsNumericalDerivative tpnd(model);

    tpnd.setParametersToDerivate(model->getUnfrozenParameters().getParameterNames());
    tpnd.enableFirstOrderDerivatives(true);
    tpnd.enableSecondOrderDerivatives(true);
    tpnd.enableSecondOrderCrossDerivatives(false); // Pseudo-Newton

    optimizer.reset(new bpp::PseudoNewtonOptimizer(&tpnd));
    optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
    optimizer->init(model->getUnfrozenParameters());
  }

  else if(options_.getOptimMethod() == "BFGS")
  {
    bpp::ThreePointsNumericalDerivative tpnd(model);

    tpnd.setParametersToDerivate(model->getUnfrozenParameters().getParameterNames());
    tpnd.enableFirstOrderDerivatives(true);
    tpnd.enableSecondOrderDerivatives(false);
    tpnd.enableSecondOrderCrossDerivatives(false);

    optimizer.reset(new bpp::BfgsMultiDimensions(&tpnd));
    optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO); // WARNING check
    optimizer->init(model->getUnfrozenParameters());
  }

  else
    throw bpp::Exception("OptimizationWrapper::Mis-specified numerical optimizer " + options_.getOptimMethod());

  bpp::FunctionStopCondition stopCond(&optimizer, options_.getFunctionTolerance());
  optimizer->setStopCondition(&stopCond);

  BackupListenerOv blo("backup_params.txt");
  optimizer->addOptimizationListener(&blo);

  std::ofstream prof("profile.txt", std::ios::out);
  std::ofstream mess("messages.txt", std::ios::out);

  bpp::StlOutputStream profiler(&prof);
  bpp::StlOutputStream messenger(&mess);

  optimizer->setProfiler(&profiler);
  optimizer->setMessageHandler(&messenger);
  
  try
    optimizer->optimize();

  catch(bpp::Exception& e)
    std::cout << "\nError during optimization, convergence might not be reached!\nmoments++ will proceed, but check log files.\n";
}

void OptimizationWrapper::computeCI() // TODO use Godambe Information Matrix
{
  // Matrix operations with BPP tools (instead of Eigen) to make our lives easier w.r.t. compatibility
  std::cout << std::endl << "Computing 95% confidence intervals of parameter estimates..." << std::endl;

  bpp::ParameterList optimParams = model -> fetchModelParameters();
  std::vector<std::string> paramNames = optimParams.getParameterNames();

  bpp::ThreePointsNumericalDerivative* tpnd = new ThreePointsNumericalDerivative(optimizedSplines);
  tpnd -> enableFirstOrderDerivatives(true);
  tpnd -> enableSecondOrderDerivatives(true);
  tpnd -> enableSecondOrderCrossDerivatives(true);

  for(size_t i = 0; i < paramNames.size(); ++i)
  {
    bpp::Constraint paramConstraint = optimParams.getParameter(paramNames[i]).getConstraint();

    // see lines 203 / 204 of ThreePointsNumericalDerivative.cpp
    double h = (1. + abs(optimParams.getParameterValue(paramNames[i]))) * tpnd -> getInterval();
    double lowerPointDerivative = optimParams.getParameterValue(paramNames[i]) - h - tpnd -> getInterval() / 2.; //conservative
    double upperPointDerivative = optimParams.getParameterValue(paramNames[i]) + h + tpnd -> getInterval() / 2.; //conservative

    if(!paramConstraint -> includes(lowerPointDerivative, upperPointDerivative))
    {
      optimParams.deleteParameter(paramNames[i]);
      std::cout << "   Numerical derivative can't be computed for " << paramNames[i] << "! Estimate is too close to boundary." << std::endl;
    }

    //updates
    paramNames = optimParams.getParameterNames();
    tpnd -> setParametersToDerivate(paramNames);

    bpp::RowMatrix< double > varCovarMatrix;
    bpp::RowMatrix< double >* hessian = bpp::NumTools::computeHessianMatrix(*tpnd, optimParams);
    bpp::MatrixTools::inv(*hessian, varCovarMatrix);

    //prints covariance matrix
    std::ofstream outMatrix("varCovar.txt", std::ios::out);
    outMatrix << paramNames[0];

    for(size_t j = 1; j < paramNames.size(); j++)
      outMatrix << " " << paramNames[j];

    outMatrix << std::endl;

    for(size_t i = 0; i < paramNames.size(); i++)
    {
      outMatrix << paramNames[i];

      for(size_t j = 0; j < paramNames.size(); j++)
        outMatrix << " " << varCovarMatrix(i, j);

      outMatrix << std::endl;
      outMatrix.close();
    }

    //print 95% CI's
    std::ofstream outCI("CI.txt", std::ios::out);

    for(size_t i = 0; i < optimParams.size(); ++i)
    {
      double upper = optimParams.getParameterValue(paramNames[i]) + 1.96 * pow(varCovarMatrix(i, i), 0.5);
      double lower = optimParams.getParameterValue(paramNames[i]) - 1.96 * pow(varCovarMatrix(i, i), 0.5);
      outCI << paramNames[i] << " " << lower << " - " << upper << std::endl;
    }

    outCI.close();
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
