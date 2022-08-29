/*
 * Author: Gustavo V. Barroso
 * Created: 29/08/2022
 * Last modified: 29/08/2022
 * Source code for moments++
 *
 */

////////////////////////////////////////////////////////////////////////////////////////////////////

#include "SumStatsLibrary.h"
#include "Mutation.h"
#include "Recombination.h"
#include "Drift.h"
#include "Selection.h"
#include "Admixture.h"
#include "SmcOptimizationWrapper.h"
#include "OptionsContainer.h"
#include "Model.h"
#include "Data.h"

int main(int argc, char *argv[]) {

  std::cout << std::endl;
  std::cout << "******************************************************************" << std::endl;
  std::cout << "*                moments++, version 0.0.1                        *" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "*            Recombination                                       *" << std::endl;
  std::cout << "*            A mosaic in the genome                              *" << std::endl;
  std::cout << "*            Endless ancestors                                   *" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "* Authors: G. Barroso                    Last Modif. 29/Oct/2022 *" << std::endl;
  std::cout << "*          A. Ragsdale                                           *" << std::endl;
  std::cout << "******************************************************************" << std::endl;
  std::cout << std::endl;

  if(argc == 1)
  {
    std::cout << "To use moments++, please fill in the params file and simply call it from the command line: moments++ params=[params_file].bpp" << std::endl;
    std::cout << "For more information, please email gvbarroso@gmail.com " << std::endl;
    return(0);
  }

  //Reads params file
  BppApplication momentspp(argc, argv, "moments++");
  momentspp.startTimer();
  map< string, string > params = momentspp.getParams();
  shared_ptr< OptionsContainer > smcOptions = make_shared< OptionsContainer >(params);

  NUMBER_OF_AVAILABLE_THREADS = smcOptions -> getNumberOfThreads();

  //Handles input sequence data
  shared_ptr< PolymorphismData > dataSet = make_shared< PolymorphismData >(smcOptions);

  std::cout << "Setting up optimisation." << std::endl;

    unsigned int numIntervals = smcOptions -> getNumberOfIntervals();
    hmmLib = make_shared< HmmStatesLibrary >(numIntervals, paramScalings, parameterAlphabet);

    //Builds MultiplePsmc object for optimisation:
    auto psmcVectorOptim = vector< shared_ptr< MmPsmc > >();
    psmcVectorOptim.reserve(dataSet -> getNumberOfDiploids());

    model = make_shared< MarkovModulatedSmc >(smcOptions, dataSet -> getSnpCalling(), hmmLib, paramScalings, categoryTransitions, parameterAlphabet);

    transitions = make_shared< MmSmcTransitionProbabilities >(model);
    emissions = make_shared< MmSmcEmissionProbabilities >(model, dataSet -> getSnpCalling());

    for(size_t i = 0; i < dataSet -> getNumberOfDiploids(); ++i) {

      shared_ptr< MmPsmc > biHaploid = make_shared< MmPsmc >(dataSet -> getNames()[i],
                                                             dataSet -> getSnpCalling()[i],
                                                             dataSet -> getBreakpoints(),
                                                             smcOptions -> getMissingBlocksLength(),
                                                             model, transitions, emissions);
      psmcVectorOptim.push_back(biHaploid);
    }

    multiMmPsmc = make_shared< MultipleMmPsmc >(psmcVectorOptim);

    SmcOptimizationWrapper smcWrapper(model, emissions, transitions, multiMmPsmc, smcOptions);

    //decides whether or not to look for backup parameters in the working dir
    if(smcOptions -> resumeOptim()) {
      ParameterList backupParams = MarkovModulatedSmc::readParametersFromFile(smcOptions -> getLabel() + "_backup_params.txt");
      smcWrapper.optimizeParameters(backupParams);
    }
    else {
      smcWrapper.optimizeParameters();
    }

    if(smcOptions -> computeCI()) {

      std::cout << std::endl << "Computing 95% confidence intervals of parameter estimates..." << std::endl;

      SplinesModel* optimizedSplines = smcWrapper.selectBestModel().get();
      ParameterList optimParams = optimizedSplines -> fetchModelParameters();
      vector< string > paramNames = optimParams.getParameterNames();

      ThreePointsNumericalDerivative* tpnd = new ThreePointsNumericalDerivative(optimizedSplines);
      tpnd -> enableFirstOrderDerivatives(true);
      tpnd -> enableSecondOrderDerivatives(true);

      if(smcOptions -> computeCovar()) {
        tpnd -> enableSecondOrderCrossDerivatives(true);
      }

      else {
        tpnd -> enableSecondOrderCrossDerivatives(false);
      }

      for(size_t i = 0; i < paramNames.size(); ++i) { //loops over paramNames because optimParams potentially changes size

        shared_ptr<Constraint> paramConstraint = optimParams.getParameter(paramNames[i]).getConstraint();

        //cf lines 203 / 204 of ThreePointsNumericalDerivative.cpp
        double h = (1. + abs(optimParams.getParameterValue(paramNames[i]))) * tpnd -> getInterval();
        double lowerPointDerivative = optimParams.getParameterValue(paramNames[i]) - h - tpnd -> getInterval() / 2.; //conservative
        double upperPointDerivative = optimParams.getParameterValue(paramNames[i]) + h + tpnd -> getInterval() / 2.; //conservative

        if(!paramConstraint -> includes(lowerPointDerivative, upperPointDerivative)) {
          optimParams.deleteParameter(paramNames[i]);
          std::cout << "   Numerical derivative can't be computed for " << paramNames[i] << "! Estimate is too close to boundary." << std::endl;
        }
      }

      //updates
      paramNames = optimParams.getParameterNames();
      tpnd -> setParametersToDerivate(paramNames);

      RowMatrix< double > varCovarMatrix;

      if(smcOptions -> computeCovar()) {

        RowMatrix< double >* hessian = NumTools::computeHessianMatrix(*tpnd, optimParams);
        MatrixTools::inv(*hessian, varCovarMatrix);
        //prints covariance matrix
        ofstream outMatrix("varCovarMatrix.txt", ios::out);
        outMatrix << paramNames[0];

        for(size_t j = 1; j < paramNames.size(); j++) {
          outMatrix << " " << paramNames[j];
        }
        outMatrix << std::endl;

        for(size_t i = 0; i < paramNames.size(); i++) {

          outMatrix << paramNames[i];

          for(size_t j = 0; j < paramNames.size(); j++) {

            outMatrix << " " << varCovarMatrix(i, j);
          }
          outMatrix << std::endl;
        }

        outMatrix.close();
      }

      else { //else the variance-covariance matrix is diagonal (ie only variances are computed)

        RowMatrix< double >* hessian = computeDiagonalHessian(*tpnd, optimParams);
        invertDiagonalMatrix(*hessian, varCovarMatrix);
      }

      //print 95% CI's
      ofstream outCI("ConfidenceIntervals.txt", ios::out);

      for(size_t i = 0; i < optimParams.size(); ++i) {

        double upper = optimParams.getParameterValue(paramNames[i]) + 1.96 * pow(varCovarMatrix(i, i), 0.5);
        double lower = optimParams.getParameterValue(paramNames[i]) - 1.96 * pow(varCovarMatrix(i, i), 0.5);
        outCI << paramNames[i] << " " << lower << " - " << upper << std::endl;
      }

      outCI.close();
    }

  momentspp.done();

  } catch(exception& e) {
    std::cout << "moments++ terminated because of an error." << std::endl;
    std::cout << e.what() << std::endl;
    return 1;
  }

  return 0;
}
