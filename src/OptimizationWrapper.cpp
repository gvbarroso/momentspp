/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 06/09/2022
 *
 */


#include <iostream>
#include <memory>

#include "SumStatsLibrary.hpp"
#include "OptimizationWrapper.hpp"
#include "BackupListenerOv.hpp"
#include "Drift.hpp"
#include "Mutation.hpp"
#include "Recombination.hpp"
#include "Migration.hpp"
#include "Epoch.hpp"
#include "Model.hpp"


void OptimizationWrapper::optimize(const SumStatsLibrary& sslib)
{
  size_t numPops = sslib.getNumPops();
  size_t minNumEpochs = options_.getMinNumberOfEpochs();
  size_t maxNumEpochs = options_.getMaxNumberOfEpochs();

  for(size_t i = minNumEpochs; i <= maxNumEpochs; ++i) // one model per count of epochs
  {
    std::string name = bpp::ToString(i) + "_epochs_model";

    std::vector<std::shared_ptr<Epoch>> epochs(i);

    bpp::ParameterList mutPl;
    bpp::ParameterList recPl;

    // the range of values that our "small" rates are allowed to take in
    std::shared_ptr<bpp::IntervalConstraint> ic = std::make_shared<bpp::IntervalConstraint>(0., 1e-5, true, true);

    mutPl.addParameter(new bpp::Parameter("r_0", 1e-8, ic));
    recPl.addParameter(new bpp::Parameter("u_0", 1e-8, ic));

    std::shared_ptr<Recombination> recOp = std::make_shared<Recombination>(rec, sslib);
    std::shared_ptr<Mutation> mutOp = std::make_shared<Mutation>(mut, sslib);

    for(size_t j = 0; j < i; ++j) // for each epoch, from past to present
    {
      std::string id = "e_" + bpp::ToString(j); // for setting the namespace for params within each epoch

      bpp::ParameterList driftPl;
      bpp::ParameterList migPl;

      for(size_t k = 0; k < numPops; ++k) // NOTE could be numPops[j]...
      {
        driftPl.addParameter(new bpp::Parameter("N_" + bpp::TextTools::ToString(k), 1e+4, bpp::Parameter::R_PLUS_STAR));

        for(size_t l = 0; l < numPops; ++l)
          migPl.addParameter(new bpp::Parameter("m_" + bpp::TextTools::ToString(k) + bpp::TextTools::ToString(l), 1e-8, ic));
      }

      std::shared_ptr<Drift> driftOp = std::make_shared<Drift>(driftPl, sslib);
      std::shared_ptr<Migration> migOp = std::shared_ptr<Migration>(migPl, sslib);

      std::vector<std::shared_ptr<Operator>> operators(4);

      // include operators in the correct order for matrix operations
      operators[0] = migOp;
      operators[1] = driftOp;
      operators[2] = recOp;
      operators[3] = mutOp;

      // define start and end of epochs as quantiles of the exp dist?
      size_t start = j * (totalNumberOfGenerations / i);// in units of generations
      size_t end = (j + 1) * (totalNumberOfGenerations / i); // in units of generations

      epochs.emplace_back(std::make_shared<Epoch>(operators, start, end, id));
    }

    Model* model = new Model(name, epochs, sslib);

    /*
     * decides whether to freeze parameters
     * eg, model->freezeParamter("r_0")
     */

    fitModel_(model);
    writeEstimatesToFile_(model);

    delete model;
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

  bpp::FunctionStopCondition stopCond(optimizer.get(), options_.getFunctionTolerance());
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

  file << "CLL = " << model->comLogLikelihood() << std::endl << std::endl;
  model->getParameterList().printParameters(file);

  file.close();
}
