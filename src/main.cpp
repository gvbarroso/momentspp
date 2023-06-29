/*
 * Author: Gustavo V. Barroso
 * Created: 29/08/2022
 * Last modified: 14/06/2023
 * Source code for moments++
 *
 */

////////////////////////////////////////////////////////////////////////////////////////////////////

#include "SumStatsLibrary.hpp"
#include "Mutation.hpp"
#include "Recombination.hpp"
#include "Drift.hpp"
#include "Selection.hpp"
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
  std::cout << "*                 moments++  version 0.0.1                       *" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "*                   \"Barrilete Cosmico\"                          *" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "*            Two-site recursions                                 *" << std::endl;
  std::cout << "*            Unraveling history                                  *" << std::endl;
  std::cout << "*            Moment by moment                                    *" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "* Authors: G. Barroso                    Last Modif. 29/Jun/2023 *" << std::endl;
  std::cout << "*          A. Ragsdale                                           *" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "******************************************************************" << std::endl;

  std::cout << "\nCompiled on: " << __DATE__ << std::endl;
  std::cout << "Compiled at: " << __TIME__ << std::endl << std::endl;

  if(argc == 1)
  {
    std::cout << "To use moments++, fill in a text file with the following options and execute from the command line:\nmomentspp params=file_name\n\n";

    std::cout << "label = # optional string, default = 'moments++'\n";
    std::cout << "demes_file = # mandatory, relative path to file in Demes format that specifies the (starting) model\n";
    std::cout << "stats_file = # optional, relative path to file listing observed summary statistics from sampled populations, default = 'none'\n\n";

    std::cout << "optimizer = # optional string, default = 'BFGS'\n";
    std::cout << "tolerance = # optional double, default = 1e-6\n";
    std::cout << "compress_moments = # optional boolean, default = TRUE\n";
    std::cout << "num_threads = # optional unsigned int, default = num_cores / 2\n";

    std::cout << "For more information, please email gvbarroso@gmail.com " << std::endl;
    return(0);
  }

  bpp::BppApplication momentspp(argc, argv, "moments++");
  momentspp.startTimer();
  std::map<std::string, std::string> params = momentspp.getParams();

  OptionsContainer options(params);

  std::cout << "\nmoments++ is using " << options.getNumThreads() << " threads.\n";
  Eigen::setNbThreads(options.getNumThreads());

  Demes demes(options.getDemesFilePath());

  size_t numEpochs = demes.getNumEpochs();
  std::vector<std::shared_ptr<Epoch>> epochs(0);
  epochs.reserve(numEpochs);

  for(size_t i = 0; i < numEpochs; ++i) // for each epoch, from past to present
  {
    std::string id = "e_" + bpp::TextTools::toString(i);

    size_t start = demes.getPopsVec()[i].front()->getStartTime(); // shared by all pops in epoch i
    size_t end = demes.getPopsVec()[i].front()->getEndTime(); // shared by all pops in epoch i

    SumStatsLibrary sslib(demes.getPopsVec()[i], options.getFactorOrder(), options.compressMoments());

    std::vector<std::shared_ptr<AbstractOperator>> operators(0);

    /* Epoch-specific operators (concern populations present in each epoch, hence parameters must follow suit)
     * Must have epoch-specific recombination and mutation operators because they depend on pop indices,
     * even though inside Model we alias r and mu across epochs
     * NOTE Admixture is modeled as the only operator in an epoch of 1 generation
     */

    if((start - end) == 1)
    {
      throw bpp::Exception("Attempted to buid Admixture operator under selection model!");

      /*if(!demes.getPulse(i).isZero(0))
      {
        operators.push_back(std::make_shared<Admixture>(demes.getPulse(i), sslib));
        //operators.back()->printTransitionLDMat(options.getLabel() + "_" + id + "_admix.csv", sslib);
      }

      else
        throw bpp::Exception("Zero Admixture matrix assigned to 1-generation Epoch!");*/
    }

    else
    {
      if(demes.getPulse(i).isZero(0))
      {
        std::shared_ptr<bpp::IntervalConstraint> ic = std::make_shared<bpp::IntervalConstraint>(0., 1e-2, true, true);
        std::shared_ptr<bpp::IntervalConstraint> icSel = std::make_shared<bpp::IntervalConstraint>(-1e-2, 0., true, true);

        std::vector<double> drift(0);
        drift.reserve(demes.getPopsVec()[i].size());

        // from (diploid) population sizes (N_j, not 2N_j) to drift parameters
        for(size_t j = 0; j < demes.getPopsVec()[i].size(); ++j)
          drift.emplace_back(1. / (2. * demes.getPopsVec()[i][j]->getSize()));

        std::shared_ptr<Drift> driftOp = std::make_shared<Drift>(drift, ic, sslib);
        std::shared_ptr<Recombination> recOp = std::make_shared<Recombination>(demes.getRec(i), ic, sslib);
        std::shared_ptr<Mutation> mutOp = std::make_shared<Mutation>(demes.getMu(i), ic, sslib);
        std::shared_ptr<Selection> selOp = std::make_shared<Selection>(demes.getSelCoeff(i), icSel, sslib);

        // only *allow* model to optimize mig params in epochs where the demes model has non-zero mig
        /*if((demes.getNumPops(i) > 1) && (!demes.getMig(i).isZero()))
        {
          operators.push_back(std::make_shared<Migration>(demes.getMig(i), ic, sslib));
          operators.back()->printDeltaLDMat(options.getLabel() + "_" + id + "_mig.csv", sslib);
        }*/

        operators.push_back(driftOp);
        operators.push_back(recOp);
        operators.push_back(mutOp);
        operators.push_back(selOp);

        for(size_t x = 0; x < operators.size(); ++x)
          operators[x]->printDeltaLDMat(options.getLabel() + "_" + id + "_op_" + bpp::TextTools::toString(x) + ".csv");

        // if previous epoch is an Admixture epoch, we correct for the 1-gen by incrementing start
        if(epochs.size() > 1 && epochs.back()->duration() == 1)
          ++start;
      }

      else
        throw bpp::Exception("Non-Zero Admixture matrix assigned to multi-generation Epoch!");
    }

    epochs.emplace_back(std::make_shared<Epoch>(id, sslib, start, end, operators, demes.getPopsVec()[i]));
    epochs.back()->printRecursions(std::cout);
    epochs.back()->printTransitionMat(options.getLabel() + "_" + id + "_transitions.csv");
  }

  epochs.front()->computeSteadyState(); // only need to have steady state in the deepest epoch

  std::ofstream fs(options.getLabel() + "_steady_state.txt");
  fs << epochs.front()->getSteadyState() << "\n";
  fs.close();

  try
  {
    if(options.getDataFilePath() == "none")
    {
      std::cout << "\nNo stats_file provided, moments++ will output expectations for input parameters.\n";
      std::shared_ptr<Model> model = std::make_shared<Model>(options.getLabel(), epochs);

      //model->getParameters().printParameters(std::cout);
      model->getIndependentParameters().printParameters(std::cout);
      model->computeExpectedSumStats();

      std::string file = model->getName() + "_expectations.txt";
      std::ofstream fout(file);

      model->printAliasedMoments(fout);

      fout.close();
      std::cout << "Check " << file << ".\n\n";
    }

    else
    {
      std::cout << "\nstats_file provided, moments++ will optimize parameters for input data.\n";
      std::shared_ptr<Data> data = std::make_shared<Data>(options.getDataFilePath());
      std::shared_ptr<Model> model = std::make_shared<Model>(options.getLabel(), epochs, data);

      OptimizationWrapper optimizer(options);
      optimizer.fitModel(model.get());
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
