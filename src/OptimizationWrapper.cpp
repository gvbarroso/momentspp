/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 01/04/2025
 *
 */

#include <iostream>
#include <memory>

#include <Bpp/Text/TextTools.h>

#include "SumStatsLibrary.hpp"
#include "OptimizationWrapper.hpp"
#include "Drift.hpp"
#include "Mutation.hpp"
#include "Recombination.hpp"
#include "Migration.hpp"
#include "Epoch.hpp"
#include "Model.hpp"

void OptimizationWrapper::fitModel(std::shared_ptr<Model> model)
{
  std::cout << "\nOptimizing model with the following parameters:\n";
  model->getUnfrozenParameters().printParameters(std::cout);

  std::unique_ptr<bpp::OptimizerInterface> optimizer;

  if(options_.getOptimMethod() == "Powell")
  {
    std::shared_ptr<bpp::ReparametrizationFunctionWrapper> rfw =
    std::make_shared<bpp::ReparametrizationFunctionWrapper>(model, model->getUnfrozenParameters());

    optimizer.reset(new bpp::PowellMultiDimensions(rfw));
    optimizer->init(rfw->getParameters());
  }

  else if(options_.getOptimMethod() == "BFGS")
  {
    std::shared_ptr<bpp::ThreePointsNumericalDerivative> tpnd =
    std::make_shared<bpp::ThreePointsNumericalDerivative>(model);

    tpnd->setParametersToDerivate(model->getUnfrozenParameters().getParameterNames());
    tpnd->enableFirstOrderDerivatives(true);
    tpnd->enableSecondOrderDerivatives(false);
    tpnd->enableSecondOrderCrossDerivatives(false);

    optimizer.reset(new bpp::BfgsMultiDimensions(tpnd));
    optimizer->setConstraintPolicy(bpp::AutoParameter::CONSTRAINTS_AUTO);
    optimizer->init(model->getUnfrozenParameters());
  }

  else
    throw bpp::Exception("OptimizationWrapper::Mis-specified numerical optimizer " + options_.getOptimMethod());

  /* NOTE this needs to be ported to bpp-core 3.0.0
  std::unique_ptr<std::ofstream> prof = std::make_unique<std::ofstream>("profile.txt", std::ios::out);
  std::unique_ptr<std::ofstream> mess = std::make_unique<std::ofstream>("messages.txt", std::ios::out);

  std::unique_ptr<std::ostream> prof_ptr = std::move(prof);
  std::unique_ptr<std::ostream> mess_ptr = std::move(mess);

  std::shared_ptr<bpp::StlOutputStream> profiler = std::make_shared<bpp::StlOutputStream>(prof_ptr);
  std::shared_ptr<bpp::StlOutputStream> messenger = std::make_shared<bpp::StlOutputStream>(mess_ptr);

  optimizer->setProfiler(profiler);
  optimizer->setMessageHandler(messenger);

  std::shared_ptr<bpp::OptimizationStopCondition> stopCond;// = std::make_shared<bpp::OptimizationStopCondition>();
  stopCond->setOptimizer(optimizer.get());
  stopCond->setTolerance(options_.getTolerance());
  optimizer->setStopCondition(stopCond);
  */

  try
  {
    optimizer->optimize();
  }

  catch(bpp::Exception& e)
  {
    std::cout << "\nError during optimization, convergence might not be reached!\nmoments++ will proceed, but check log files.\n";
  }
}

/*void OptimizationWrapper::computeCI()
{
  // TODO use Godambe Information Matrix
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
}*/

void OptimizationWrapper::writeEstimatesToFile_(std::shared_ptr<Model> model)
{
  std::ofstream file;
  file.open(model->getName() + "_estimates.txt");

  file << "CLL = " << - model->getValue() << "\n\n";

  for(size_t i = 0; i < model->getEpochs().size(); ++i)
  {
    std::shared_ptr<Epoch> epoch = model->getEpochs()[i];
    file << epoch->getNamespace() << "\t" << epoch->start() << "\t" << epoch->end() << "\n";

    for(size_t j = 0; j < epoch->getParameters().size(); ++j)
      file << epoch->getParameters()[i].getName() << "=" << epoch->getParameters()[i].getValue() << "\n";

    file << "\n";
  }

  file.close();
}
