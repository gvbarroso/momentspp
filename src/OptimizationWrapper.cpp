/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 15/02/2023
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


void OptimizationWrapper::optimize(const Data& data, const Demes& demes)
{
  size_t numEpochs = demes.getNumEpochs();
  std::string modelName = options_.getLabel();

  std::vector<std::shared_ptr<Epoch>> epochs(0);
  epochs.reserve(numEpochs);

  // the range of values that our "small" rates are allowed to take in
  std::shared_ptr<bpp::IntervalConstraint> ic = std::make_shared<bpp::IntervalConstraint>(0., 1e-3, true, true);

  for(size_t i = 0; i < numEpochs; ++i) // for each epoch, from past to present
  {
    std::string id = "e_" + bpp::TextTools::toString(i); // for setting the namespace for params within each epoch

    size_t start = i * (options_.getTotalNumberOfGenerations() / numEpochs);
    size_t end = (i + 1) * (options_.getTotalNumberOfGenerations() / numEpochs);

    SumStatsLibrary sslib(options_.getOrder(), demes.getPopMaps()[i]);

    /* Epoch-specific operators (concern populations present in each epoch, hence parameters must follow suit)
     * must have epoch-specific recombination and mutation operators because they depend on pop indices (popMaps[i]),
     * even though we prob. want single r and mu params in Model --> alias r and mu across epochs?
     */

    std::vector<double> drift = options_.getInitPopSizes(); // should get this from Demes instead
    for(size_t j = 0; j < drift.size(); ++j) // from (diploid) population sizes (N_j, not 2N_j) to (diploid) drift parameters
      drift[j] = 1. / (2. * drift[j]);

    // WARNING ad-hockery for testing bottleneck
    if(i == 1)
    {
      for(size_t j = 0; j < drift.size(); ++j)
        drift[j] *= 10.;
    }

    std::shared_ptr<Drift> driftOp = std::make_shared<Drift>(drift, ic, sslib);
    std::shared_ptr<Recombination> recOp = std::make_shared<Recombination>(options_.getInitR(), ic, sslib);
    std::shared_ptr<Mutation> mutOp = std::make_shared<Mutation>(options_.getInitMu(), ic, sslib);

    std::vector<std::shared_ptr<AbstractOperator>> operators(0);
    operators.reserve(4);

    if(demes.getNumPops() > 1)
    {
      std::shared_ptr<Migration> migOp = std::make_shared<Migration>(options_.getInitMig(), ic, sslib);
      operators.emplace_back(migOp);
    }

    // include operators in the "correct" order for matrix operations
    operators.emplace_back(driftOp);
    operators.emplace_back(recOp);
    operators.emplace_back(mutOp);

    epochs.emplace_back(std::make_shared<Epoch>(id, sslib, start, end, operators, demes.getPopMaps()[i]));
  }

  Model* model = new Model(modelName, epochs, data);
  model->computeExpectedSumStats();
  model->aliasMoments();

  std::ofstream fout(modelName + "_final_unsorted.txt");
  model->printAliasedMoments(fout);
  fout.close();

  //fitModel_(model);
  //writeEstimatesToFile_(model);

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

  file << "CLL = " << model->compLogLikelihood() << "\n\n";

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
