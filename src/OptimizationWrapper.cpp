/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 13/09/2022
 *
 */

#include <iostream>
#include <memory>

#include <Bpp/Text/TextTools.h>

#include "SumStatsLibrary.hpp"
#include "OptimizationWrapper.hpp"
#include "BackupListenerOv.hpp"
#include "Drift.hpp"
#include "Mutation.hpp"
#include "Recombination.hpp"
#include "Migration.hpp"
#include "Epoch.hpp"
#include "Model.hpp"


void OptimizationWrapper::optimize()
{
  std::vector<size_t> numPops = options_.getNumbersOfPopulations(); // one per epoch bc of population splits / merges
  size_t numEpochs = options_.getNumberOfEpochs();

  std::string name = bpp::TextTools::toString(numEpochs) + "_epochs_model";

  std::vector<std::shared_ptr<Epoch>> epochs(0);
  epochs.reserve(numEpochs);

  bpp::ParameterList mutPl;
  bpp::ParameterList recPl;

  // the range of values that our "small" rates are allowed to take in
  std::shared_ptr<bpp::IntervalConstraint> ic = std::make_shared<bpp::IntervalConstraint>(0., 1e-5, true, true);

  mutPl.addParameter(new bpp::Parameter("r_0", 1e-8, ic));
  recPl.addParameter(new bpp::Parameter("u_0", 1e-8, ic));

  for(size_t i = 0; i < numEpochs; ++i) // for each epoch, from past to present
  {
    std::string id = "e_" + bpp::TextTools::toString(i); // for setting the namespace for params within each epoch

    // define start and end of epochs as quantiles of the exp dist?
    size_t start = (numEpochs - i) * (options_.getTotalNumberOfGenerations() / numEpochs); // in units of generations
    size_t end = (numEpochs - i - 1) * (options_.getTotalNumberOfGenerations() / numEpochs); // in units of generations

    bpp::ParameterList driftPl;
    bpp::ParameterList migPl;

    for(size_t j = 0; j < numPops[i]; ++j) // for each population modeled in epoch i
    {
      driftPl.addParameter(new bpp::Parameter("N_" + bpp::TextTools::toString(j), 1e+4, bpp::Parameter::R_PLUS_STAR));

      for(size_t k = 0; k < numPops[i]; ++k)
      {
        if(k != j)
          migPl.addParameter(new bpp::Parameter("m_" + bpp::TextTools::toString(j) + bpp::TextTools::toString(k), 1e-8, ic));
      }
    }

    std::map<size_t, std::pair<size_t, size_t>> popMap; // TODO read table from file, check if there's standard format in demes
    SumStatsLibrary sslib;
    sslib.initStatsVector(popMap);

    std::shared_ptr<Drift> driftOp = std::make_shared<Drift>(driftPl, sslib);
    std::shared_ptr<Migration> migOp = std::make_shared<Migration>(migPl, sslib);
    std::shared_ptr<Recombination> recOp = std::make_shared<Recombination>(recPl, sslib);
    std::shared_ptr<Mutation> mutOp = std::make_shared<Mutation>(mutPl, sslib);

    std::vector<std::shared_ptr<Operator>> operators(0);
    operators.reserve(4);

    // include operators in the correct order for matrix operations
    operators.emplace_back(migOp);
    operators.emplace_back(driftOp);
    operators.emplace_back(recOp);
    operators.emplace_back(mutOp);

    epochs.emplace_back(std::make_shared<Epoch>(operators, popMap, start, end, id));
  }

  Model* model = new Model(name, epochs, sslib);
  model->computeSteadyState();

  /*
   * decides whether to freeze parameters
   * eg, model->freezeParamter("r_0")
   */

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

  BackupListenerOv blo("backup_params.txt");
  optimizer->addOptimizationListener(&blo);

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

/*void OptimizationWrapper::computeCI() // TODO use Godambe Information Matrix
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
}*/

void OptimizationWrapper::writeEstimatesToFile_(Model* model)
{
  std::ofstream file;
  file.open(model->getName() + "_estimates.txt");

  file << "CLL = " << model->comLogLikelihood() << std::endl << std::endl;

  for(size_t i = 0; i < model->getEpochs().size(); ++i)
  {
    std::shared_patr<Epoch> epoch = model->getEpochs()[i];
    file << epoch->getPrefix() << "\t" << epoch->start() << "\t" << epoch->end() << "\t";

    for(size_t j = 0; j < epoch->getParameters().size(); ++j)
    {
      file << epoch->getParameters()[i].getParameterName() << "=" << epoch->getParameters()[i].getParameterValue();

      if(i < epoch->getParameters().size() - 1)
        file << "\t";

      else
        file << "\n";
    }
  }

  file.close();
}
