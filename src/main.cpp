/*
 * Author: Gustavo V. Barroso
 * Created: 29/08/2022
 * Last modified: 24/03/2025
 * Source code for moments++
 *
 */

////////////////////////////////////////////////////////////////////////////////////////////////////

#include "SumStatsLibrary.hpp"
#include "Mutation.hpp"
#include "Recombination.hpp"
#include "Drift.hpp"
#include "Selection.hpp"
//#include "NeutralMigration.hpp"
//#include "NeutralAdmixture.hpp"
//#include "Migration.hpp"
//#include "Admixture.hpp"
#include "OptimizationWrapper.hpp"
#include "OptionsContainer.hpp"
#include "Model.hpp"
#include "Data.hpp"
#include "Demes.hpp"

int main(int argc, char *argv[]) {

  std::cout << std::endl;
  std::cout << "******************************************************************" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "*                  moments++  version 0.1                      *" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "*            Two-site recursions                                 *" << std::endl;
  std::cout << "*            Unraveling history                                  *" << std::endl;
  std::cout << "*            Moment by moment                                    *" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "* Authors: G. V. Barroso                 Last Modif. 14/Mar/2025 *" << std::endl;
  std::cout << "*          A. P. Ragsdale                                        *" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "******************************************************************" << std::endl;

  std::cout << "\nCompiled on: " << __DATE__ << std::endl;
  std::cout << "Compiled at: " << __TIME__ << std::endl << std::endl;

  /*
   * 1. Variance in Heterozigosity across left and right loci  (p^2 * q^2)
   * 2. To compress basis by adding (averaging) rows of uncompressed Matrices, then removing corresponding row and column:
   *    https://stackoverflow.com/questions/13290395/how-to-remove-a-certain-row-or-column-while-using-eigen-library-c
   */

  if(argc == 1)
  {
    std::cout << "Usage:\n";
    std::cout << "momentspp param=opt.bpp\n\n";

    std::cout << "\nThe github repository contains instructions on how to write the options file:\n";
    std::cout << "https://github.com/gvbarroso/momentspp/tree/main/doc" << std::endl << std::endl;
    std::cout << "\nIf you have any doubts, please email gvbarroso@gmail.com " << std::endl;
    return(0);
  }

  bpp::BppApplication momentspp(argc, argv, "moments++");
  momentspp.startTimer();
  std::map<std::string, std::string> params = momentspp.getParams();

  OptionsContainer options(params);

  std::cout << "\nmoments++ is using " << options.getNumThreads() << " threads.\n";
  Eigen::setNbThreads(options.getNumThreads());

  Demes demes(options.getDemesFilePath());

  std::cout << "Assembling Operators and Epoch objects..."; std::cout.flush();

  size_t numEpochs = demes.getNumEpochs();
  std::vector<std::shared_ptr<Epoch>> epochs(0);
  epochs.reserve(numEpochs);

  std::vector<size_t> factorOrder = options.getFactorOrder(); // one value per epoch to avoid underflow e.g. after a bottleneck
  if(factorOrder.size() == 1)
  {
    for(size_t i = 1; i < numEpochs; ++i)
      factorOrder.push_back(factorOrder[0]);
  }
  else if(factorOrder.size() != numEpochs)
    throw bpp::Exception("Main::Number of Factor Orders must be either 1 or equal to the number of Epochs in the model!");

  if(std::any_of(std::begin(factorOrder), std::end(factorOrder), [](size_t x) { return x < 1; }))
    throw bpp::Exception("Main::All Factor Orders must be greather than zero!");

  auto finder = std::adjacent_find(std::begin(factorOrder), std::end(factorOrder), std::less<size_t>());
  if(finder != std::end(factorOrder))
    throw bpp::Exception("Main::Factor Orders can not increase over time!");

  for(size_t i = 0; i < numEpochs; ++i) // for each epoch, from past to present
  {
    std::string id = "e_" + bpp::TextTools::toString(i);

    size_t start = demes.getPopsVec()[i].front()->getStartTime(); // shared by all pops in epoch i
    size_t end = demes.getPopsVec()[i].front()->getEndTime(); // shared by all pops in epoch i

    SumStatsLibrary sslib(demes.getPopsVec()[i], factorOrder[i], options.compressMoments());

    std::vector<std::shared_ptr<AbstractOperator>> operators(0);

    /* Epoch-specific operators (concern populations present in each epoch, hence parameters must follow suit)
     * Must have epoch-specific recombination and mutation operators because they depend on pop indices,
     * even though inside Model class we often choose to alias r and mu across epochs and pops.
     */

    // NOTE this current implementation generates a problem if the user wants other 1-gen epochs for some reason
    if((start - end) == 1) // Admixture is modeled as the only operator in an epoch of 1 generation
    {
      /*if(!demes.getPulse(i).isZero(0))
      {
        operators.push_back(std::make_shared<Admixture>(demes.getPulse(i), sslib));
        //operators.back()->printTransitionLDMat(options.getLabel() + "_" + id + "_admix.csv", sslib);
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

        std::vector<long double> drift(0);
        drift.reserve(demes.getPopsVec()[i].size());

        // from (diploid) population sizes (N_j, not 2N_j) to drift parameters
        for(size_t j = 0; j < demes.getPopsVec()[i].size(); ++j)
          drift.emplace_back(1. / (2. * demes.getPopsVec()[i][j]->getSize()));

        std::shared_ptr<Selection> selOp = std::make_shared<Selection>(demes.getSelCoeffs(i), icSel, sslib);
        std::shared_ptr<Recombination> recOp = std::make_shared<Recombination>(demes.getRecs(i), icRec, sslib);
        std::shared_ptr<Mutation> mutOp = std::make_shared<Mutation>(demes.getLeftFactor(), demes.getMus(i), ic, sslib);
        std::shared_ptr<Drift> driftOp = std::make_shared<Drift>(drift, ic, sslib);

        /*// only *allow* model to include mig params in epochs where the demes model has non-zero mig
        if((demes.getNumPops(i) > 1) && (!demes.getMig(i).isZero()))
        {
          operators.push_back(std::make_shared<Migration>(demes.getMig(i), ic, sslib));
          //operators.back()->printDeltaLDMat(options.getLabel() + "_" + id + "_mig.csv");
        }*/

        operators.push_back(selOp);
        operators.push_back(recOp);
        operators.push_back(mutOp);
        operators.push_back(driftOp);

        if(options.verbose())
        {
          for(size_t j = 0; j < operators.size(); ++j)
            operators[j]->printDeltaLDMat(options.getLabel() + "_" + id + "_op_" + bpp::TextTools::toString(j) + ".csv");
        }

        // if immediately previous epoch is an Admixture epoch, we correct for the 1-gen by incrementing start
        if(epochs.size() > 1 && epochs.back()->duration() == 1)
          ++start;
      }

      else
        throw bpp::Exception("Main::Non-Zero Admixture matrix assigned to multi-generation Epoch!");
    }

    epochs.emplace_back(std::make_shared<Epoch>(id, sslib, start, end, operators, demes.getPopsVec()[i]));

    if(options.verbose())
    {
      epochs.back()->printRecursions(std::cout);
      epochs.back()->printTransitionMat(options.getLabel() + "_" + id + "_O_" + bpp::TextTools::toString(factorOrder[0]) +"_transitions.csv");
    }
  }

  if(options.verbose())
  {
    epochs.front()->computePseudoSteadyState();
    std::ofstream pseudo(options.getLabel() + "_O_" + bpp::TextTools::toString(factorOrder[0]) + "_pseudo_steady-state.txt");
    epochs.front()->printMoments(pseudo);
    pseudo.close();

    epochs.front()->computeEigenSteadyState();
    std::ofstream eigen(options.getLabel() + "_O_" + bpp::TextTools::toString(factorOrder[0]) + "_eigen_steady-state.txt");
    epochs.front()->printMoments(eigen);
    eigen.close();
  }

  if(options.getInitStatsFilePath() == "none")
    epochs.front()->computeEigenSteadyState(); // only need steady state in the deep-most epoch

  else
    epochs.front()->getSslib().readStatsFromFile(options.getInitStatsFilePath()); // NOTE mind Order of (1-2p) factors

  std::cout << "done.\n\nBuilding Model now.";

  try
  {
    if(options.getDataFilePath() == "none")
    {
      std::cout << "\nNo obs_stats_file provided, moments++ will\noutput expectations for input parameters.\n\n";

      std::shared_ptr<Model> model = std::make_shared<Model>(options.getLabel(), epochs);
      model->getIndependentParameters().printParameters(std::cout);
      model->computeExpectedSumStats();

      std::string fileName = model->getName() + "_O_" + bpp::TextTools::toString(factorOrder[0]) + "_expectations.txt";
      std::ofstream fout(fileName);

      model->printAliasedMoments(fout);
      if(numEpochs > 1 && options.getTimeSteps() > 0)
        model->printHetMomentsIntermediate(model->getName() + "_O_" + bpp::TextTools::toString(factorOrder[0]), options.getTimeSteps());

      fout.close();
      std::cout << "\nCheck output file " << fileName << "\n\n";
    }

    else
    {
      std::cout << "\nStats_file provided, moments++ will optimize parameters for input data.\n";

      std::shared_ptr<Data> data = std::make_shared<Data>(options.getDataFilePath());
      std::shared_ptr<Model> model = std::make_shared<Model>(options.getLabel(), epochs, data);
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

  momentspp.done();
  return 0;
}
