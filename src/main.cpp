/*
 * Author: Gustavo V. Barroso
 * Created: 29/08/2022
 * Last modified: 14/10/2025
 * Source code for moments++
 *
 */

////////////////////////////////////////////////////////////////////////////////////////////////////

#include "SumStatsLibrary.hpp"
#include "Mutation.hpp"
#include "Recombination.hpp"
#include "Drift.hpp"
#include "Selection.hpp"
// #include "NeutralMigration.hpp"
// #include "NeutralAdmixture.hpp"
// #include "Migration.hpp"
// #include "Admixture.hpp"
#include "OptimizationWrapper.hpp"
#include "OptionsContainer.hpp"
#include "Model.hpp"
#include "Data.hpp"
#include "Demes.hpp"

int main(int argc, char* argv[])
{

  std::cout << std::endl;
  std::cout << "******************************************************************" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "*                    moments++  version 0.1                      *" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "*            Two-site recursions                                 *" << std::endl;
  std::cout << "*            Unraveling history                                  *" << std::endl;
  std::cout << "*            Moment by moment                                    *" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "* Authors: G. V. Barroso                 Last Modif. 15/Oct/2025 *" << std::endl;
  std::cout << "*          A. P. Ragsdale                                        *" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "******************************************************************" << std::endl;

  std::cout << "\nCompiled on: " << __DATE__ << std::endl;
  std::cout << "Compiled at: " << __TIME__ << std::endl << std::endl;

  /*
   * NOTE there are three classes of numerical issues we need to handle:
   * 1. Multiplying "large" (O(1), Transition Matrix entries) and tiny (O(1e-12) or less, moments expectations) numbers
   * 2. Truncation of the Selection operator (moment-closure approximation) may lead to entries of TM > 1 --> worse for largest |s|
   * 3. High order of 1-2p factors combined with low population sizes lead to entries of TM < - 1 --> violation of Kingman's coalescent
   *
   * TODO homogenous/in-homogenous system (remove I moment?)
   * https://stackoverflow.com/questions/13290395/how-to-remove-a-certain-row-or-column-while-using-eigen-library-c
   *
   * TODO add method to scale matrices by 2Nanc
   *
   * 1. Variance in Heterozigosity across left and right loci  (p^2 * q^2)
   */

  if(argc == 1)
  {
    std::cout << "Usage:\n";
    std::cout << "momentspp param=opt.bpp\n\n";

    std::cout << "\nThe github repository contains instructions on how to write the options file:\n";
    std::cout << "https://github.com/gvbarroso/momentspp/tree/main/doc\n\n";
    std::cout << "\nIf you have any doubts, please email gvbarroso@gmail.com\n";
    return (0);
  }

  bpp::BppApplication momentspp(argc, argv, "moments++");
  momentspp.startTimer();
  std::map<std::string, std::string> params = momentspp.getParams();

  OptionsContainer options(params);

  if(options.highPrecision()) // if more than 16-digit precision (double)
    mpfr::mpreal::set_default_prec(mpfr::digits2bits(options.getDigits()));

  std::string scalar = "double";
  if(options.highPrecision())
    scalar = "mpfr::mpreal";

  std::cout << "\nmoments++ is using " << options.getNumThreads() << " threads.\n";
  std::cout << "numerical precision: " << options.getDigits() << " digits. Scalar type: " << scalar << "\n";

  omp_set_num_threads(options.getNumThreads());
  Eigen::setNbThreads(options.getNumThreads());
  Eigen::initParallel();

  Demes demes(options.getDemesFilePath());

  std::cout << "Assembling Operators and Epoch objects...\n\n";

  size_t numEpochs = demes.getNumEpochs();
  std::vector<std::shared_ptr<Epoch>> epochs(0);
  epochs.reserve(numEpochs);

  std::vector<size_t> factorOrder = options.getFactorOrder(); // one value per epoch to avoid underflow e.g. after a bottleneck NOTE still?
  if(factorOrder.size() == 1)
  {
    for(size_t i = 1; i < numEpochs; ++i)
      factorOrder.push_back(factorOrder[0]);
  }

  else if(factorOrder.size() != numEpochs)
    throw bpp::Exception("Main::Number of Factor Orders must be either 1 or equal to the number of Epochs in the model!");

  if(std::any_of(std::begin(factorOrder), std::end(factorOrder), [](size_t x) { return x < 1; }))
    throw bpp::Exception("Main::All Factor Orders must be greater than zero!");

  auto finder = std::adjacent_find(std::begin(factorOrder), std::end(factorOrder), std::less<size_t>());
  if(finder != std::end(factorOrder))
    throw bpp::Exception("Main::Factor Orders can not increase over time!");

  for(size_t i = 0; i < numEpochs; ++i) // for each epoch, from past to present
  {
    std::string id = "e_" + bpp::TextTools::toString(i);

    size_t start = demes.getPopsVec()[i].front()->getStartTime(); // shared by all pops in epoch i
    size_t end = demes.getPopsVec()[i].front()->getEndTime();     // shared by all pops in epoch i

    SumStatsLibrary sslib(demes.getPopsVec()[i], factorOrder[i], options.compressMoments());

    std::vector<std::shared_ptr<AbstractOperator>> operators(0);

    /* Epoch-specific operators (concern populations present in each epoch, hence parameters must
     * follow suit) Must have epoch-specific recombination and mutation operators because they
     * depend on pop indices, even though inside Model class we often choose to alias r and mu
     * across epochs and pops.
     */

    // NOTE this current implementation generates a problem if the user wants other 1-gen epochs for
    // some reason
    if((start - end) == 1) // Admixture is modeled as the only operator in an epoch of 1 generation
    {
      /*if(!demes.getPulse(i).isZero(0))
      {
        operators.push_back(std::make_shared<Admixture>(demes.getPulse(i), sslib, highPrec));
        //operators.back()->printTransitionLDMat(options.getLabel() + "_" + id + "_admix.csv",
      sslib);
      }

      else
        throw bpp::Exception("Main::Zero Admixture matrix assigned to 1-generation Epoch!");*/
    }

    else
    {
      if(demes.getPulse(i).isZero(0))
      {
        std::shared_ptr<bpp::IntervalConstraint> ic = std::make_shared<bpp::IntervalConstraint>(0., 1e-2, true, true);
        std::shared_ptr<bpp::IntervalConstraint> icRec = std::make_shared<bpp::IntervalConstraint>(0., 0.5 + 1e-6, true, true);
        std::shared_ptr<bpp::IntervalConstraint> icSel = std::make_shared<bpp::IntervalConstraint>(-1e-2, 0., true, true);

        std::vector<double> drift(0);
        drift.reserve(demes.getPopsVec()[i].size());

        // from (diploid) population sizes (N_j, not 2N_j) to drift parameters
        for(size_t j = 0; j < demes.getPopsVec()[i].size(); ++j)
          drift.emplace_back(1. / (2. * demes.getPopsVec()[i][j]->getSize()));

        bool highPrec = options.highPrecision(); // is high precision?
        std::shared_ptr<Selection> selOp = std::make_shared<Selection>(demes.getSelCoeffs(i), icSel, sslib, highPrec);
        std::shared_ptr<Recombination> recOp = std::make_shared<Recombination>(demes.getRecs(i), icRec, sslib, highPrec);
        std::shared_ptr<Mutation> mutOp = std::make_shared<Mutation>(demes.getLeftFactor(), demes.getMus(i), ic, sslib, highPrec);
        std::shared_ptr<Drift> driftOp = std::make_shared<Drift>(drift, ic, sslib, highPrec);

        /*// only *allows* model to include mig params in epochs where the demes model has non-zero
        mig if((demes.getNumPops(i) > 1) && (!demes.getMig(i).isZero()))
        {
          operators.push_back(std::make_shared<Migration>(demes.getMig(i), ic, sslib, highPrec));
          //operators.back()->printDeltaLDMat(options.getLabel() + "_" + id + "_mig.csv");
        }*/

        operators.push_back(selOp);
        operators.push_back(recOp);
        operators.push_back(mutOp);
        operators.push_back(driftOp);

        if(options.verbose()) // logs "delta" matrices for each operator
        {
          for(size_t j = 0; j < operators.size(); ++j)
            operators[j]->printDeltaLDMat(options.getLabel() + "_" + id + "_O_" +
                                          bpp::TextTools::toString(factorOrder[0]) + "_op_" +
                                          bpp::TextTools::toString(j) + ".csv");
        }

        // if immediately previous epoch is an Admixture epoch, we correct for the 1-gen by incrementing start
        if(epochs.size() > 1 && epochs.back()->duration() == 1) // NOTE epochs.back() not initialized yet
          ++start;
      }

      else
        throw bpp::Exception("Main::Non-Zero Admixture matrix assigned to multi-generation Epoch!");
    }

    // time flows from left to right, with epoch[0] (epoch.front()) => most ancient epoch
    epochs.emplace_back(std::make_shared<Epoch>(id, sslib, start, end, options.multiplyOperators(), demes.getPopsVec()[i], operators));

    if(options.verbose()) // logs transition matrices and steady-states for each epoch
    {
      epochs.back()->printRecursions(std::cout);
      epochs.back()->printTransitionMat(options.getLabel() + "_" + id + "_O_" + bpp::TextTools::toString(factorOrder[0]) + "_transitions.csv");
      epochs.back()->printConditionNumber();

      epochs.back()->computeEigenSteadyState();
      std::ofstream eigen(options.getLabel() + "_" + id + "_O_" + bpp::TextTools::toString(factorOrder[0]) + "_eigen_steady-state.txt");
      epochs.back()->printMoments(eigen);
      eigen.close();

      if(options.continuousTime())
      {
        epochs.back()->computePowerSteadyStateContinuous();
        std::ofstream power(options.getLabel() + "_" + id + "_O_" + bpp::TextTools::toString(factorOrder[0]) + "_power_steady-state.txt");
        epochs.back()->printMoments(power);
        power.close();
      }

      else
      {
        epochs.back()->computePowerSteadyStateDiscrete();
        std::ofstream power(options.getLabel() + "_" + id + "_O_" + bpp::TextTools::toString(factorOrder[0]) + "_power_steady-state.txt");
        epochs.back()->printMoments(power);
        power.close();
      }
    }
  } // ends loop over epochs

  // WARNING: this if clause is for speed in testing TODO remove it
  if(!options.verbose())
  {
  std::cout << "\nDone with Operators.\nComputing steady state..."; std::cout.flush();

  if(options.getInitStatsFilePath() == "none") // only need steady state in the deep-most epoch (epoch.front())
  {
    if(options.getSteadyStateMethod() == "eigen")
      epochs.front()->computeEigenSteadyState();

    else if(options.getSteadyStateMethod() == "power")
    {
      if(options.continuousTime())
        epochs.front()->computePowerSteadyStateContinuous();

      else
        epochs.front()->computePowerSteadyStateDiscrete();
    }

    else
      throw bpp::Exception("Main::Mis-specified steady-state method (should be 'eigen' or 'power': " + options.getSteadyStateMethod());
  }

  else
    epochs.front()->getSslib().readStatsFromFile(options.getInitStatsFilePath()); // NOTE mind Order of (1-2p) factors

  std::cout << "done.\n\nBuilding Model now.";

  try
  {
    if(options.getDataFilePath() == "none")
    {
      std::cout << "\nNo obs_stats_file provided, moments++ will\noutput expectations for input parameters.\n\n";

      std::shared_ptr<Model> model = std::make_shared<Model>(options.getLabel(),
                                                             epochs,
                                                             options.continuousTime(),
                                                             options.getDt(),
                                                             options.getTotalTimeIntegration(),
                                                             options.getToleranceIntegration());
      model->getIndependentParameters().printParameters(std::cout);
      model->computeExpectedSumStats();

      std::string fileName = model->getName() + "_O_" + bpp::TextTools::toString(factorOrder[0]) + "_expectations.txt";
      std::ofstream fout(fileName);
      model->printAliasedMoments(fout);
      fout.close();

      if(numEpochs > 1 && options.getTimeSteps() > 0)
        model->printMomentsIntermediate(model->getName() + "_O_" + bpp::TextTools::toString(factorOrder[0]),
                                        options.getTimeSteps(),
                                        options.getMomNamesIntermediate());

      std::cout << "\nCheck output file " << fileName << "\n\n";
    }

    else
    {
      std::cout << "\nStats_file provided, moments++ will optimize parameters for input data.\n";

      std::shared_ptr<Data> data = std::make_shared<Data>(options.getDataFilePath());
      std::shared_ptr<Model> model = std::make_shared<Model>(options.getLabel(),
                                                             epochs,
                                                             data,
                                                             options.continuousTime(),
                                                             options.getDt(),
                                                             options.getTotalTimeIntegration(),
                                                             options.getToleranceIntegration());

      model->compressParameters(options.aliasEpochsParams(), options.aliasPopsParams());

      std::cout << "\n\nList of parameters to be optimized:\n";
      model->getIndependentParameters().printParameters(std::cout);

      OptimizationWrapper optimizer(options);
      optimizer.fitModel(model);
    }
  }

  catch(std::exception& e)
  {
    std::cout << "moments++ terminated because of an error!" << std::endl;
    std::cout << e.what() << std::endl;

    return 1;
  }
  } // NOTE end if clause for speed in testing

  momentspp.done();
  return 0;
}
