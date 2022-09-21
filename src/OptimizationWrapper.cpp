/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 21/09/2022
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


void OptimizationWrapper::optimize(const PolymorphismData& data)
{
  std::vector<size_t> popIndices(0);
  size_t numEpochs = popList_.size();

  std::string name = bpp::TextTools::toString(numEpochs) + "_epochs_model";

  std::vector<std::shared_ptr<Epoch>> epochs(0);
  epochs.reserve(numEpochs);

  // the range of values that our "small" rates are allowed to take in
  std::shared_ptr<bpp::IntervalConstraint> ic = std::make_shared<bpp::IntervalConstraint>(0., 1e-5, true, true);

  // for now, all epochs share recombination and mutation parameters
  for(size_t i = 0; i < numEpochs; ++i) // for each epoch, from past to present
  {
    std::string id = "e_" + bpp::TextTools::toString(i); // for setting the namespace for params within each epoch

    // define start and end of epochs as quantiles of the exp dist?
    size_t start = (numEpochs - i) * (options_.getTotalNumberOfGenerations() / numEpochs); // in units of generations
    size_t end = (numEpochs - i - 1) * (options_.getTotalNumberOfGenerations() / numEpochs); // in units of generations

    // parse populations file
    std::map<size_t, std::shared_ptr<Population>> pops;
    SumStatsLibrary sslib(2, popList_[i]);

    // Epoch-specific (w.r.t populations present, hence parameters) operators
    std::shared_ptr<Drift> driftOp = std::make_shared<Drift>(sslib);
    std::shared_ptr<Migration> migOp = std::make_shared<Migration>(ic, sslib);
    // must have epoch-specific recombination and mutation operators because they depend on pop indices (popMaps[i]),
    // even though we prob. want single r and mu params in Model
    std::shared_ptr<Recombination> recOp = std::make_shared<Recombination>(ic, sslib);
    std::shared_ptr<Mutation> mutOp = std::make_shared<Mutation>(ic, sslib);

     // include operators in the correct order for matrix operations
    std::vector<std::shared_ptr<AbstractOperator>> operators(0);
    operators.reserve(4);
    operators.emplace_back(migOp);
    operators.emplace_back(driftOp);
    operators.emplace_back(recOp);
    operators.emplace_back(mutOp);

    epochs.emplace_back(std::make_shared<Epoch>(sslib, start, end, id, operators, pops));
  }

  Model* model = new Model(name, epochs, data);
  // TODO alias r and mu among epochs

  fitModel_(model);
  writeEstimatesToFile_(model);

  delete model;
}

void OptimizationWrapper::fitModel_(Model* model)
{
  std::cout << "\nOptimizing model with the following parameters:\n";
  model->getUnfrozenParameters().printParameters(std::cout);

  std::unique_ptr<bpp::Optimizer> optimizer;

  if(options_.getOptimMethod() == "Powell")
  {
    bpp::ReparametrizationFunctionWrapper rfw(model, model->getUnfrozenParameters());

    optimizer.reset(new bpp::PowellMultiDimensions(&rfw));
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
    optimizer->setConstraintPolicy(bpp::AutoParameter::CONSTRAINTS_AUTO);
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
    optimizer->setConstraintPolicy(bpp::AutoParameter::CONSTRAINTS_AUTO);
    optimizer->init(model->getUnfrozenParameters());
  }

  else
    throw bpp::Exception("OptimizationWrapper::Mis-specified numerical optimizer " + options_.getOptimMethod());

  std::unique_ptr<bpp::FunctionStopCondition> stopCond = std::make_unique<bpp::FunctionStopCondition>(optimizer.get(), options_.getTolerance());
  optimizer->setStopCondition(*stopCond);

  std::ofstream prof("profile.txt", std::ios::out);
  std::ofstream mess("messages.txt", std::ios::out);

  bpp::StlOutputStream profiler(&prof);
  bpp::StlOutputStream messenger(&mess);

  optimizer->setProfiler(&profiler);
  optimizer->setMessageHandler(&messenger);
  
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

void OptimizationWrapper::writeEstimatesToFile_(Model* model)
{
  std::ofstream file;
  file.open(model->getName() + "_estimates.txt");

  file << "CLL = " << model->comLogLikelihood() << std::endl << std::endl;

  for(size_t i = 0; i < model->getEpochs().size(); ++i)
  {
    std::shared_ptr<Epoch> epoch = model->getEpochs()[i];
    file << epoch->getNamespace() << "\t" << epoch->start() << "\t" << epoch->end() << "\t";

    for(size_t j = 0; j < epoch->getParameters().size(); ++j)
    {
      file << epoch->getParameters()[i].getName() << "=" << epoch->getParameters()[i].getValue();

      if(i < epoch->getParameters().size() - 1)
        file << "\t";

      else
        file << "\n";
    }
  }

  file.close();
}
