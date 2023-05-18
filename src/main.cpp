/*
 * Author: Gustavo V. Barroso
 * Created: 29/08/2022
 * Last modified: 02/05/2023
 * Source code for moments++
 *
 */

////////////////////////////////////////////////////////////////////////////////////////////////////

#include "SumStatsLibrary.hpp"
#include "Mutation.hpp"
#include "Recombination.hpp"
#include "Drift.hpp"
#include "Migration.hpp"
//#include "Selection.hpp"
#include "Admixture.hpp"
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
  std::cout << "* Authors: G. Barroso                    Last Modif. 18/May/2023 *" << std::endl;
  std::cout << "*          A. Ragsdale                                           *" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "******************************************************************" << std::endl;

  std::cout << "\nCompiled on: " << __DATE__ << std::endl;
  std::cout << "Compiled at: " << __TIME__ << std::endl << std::endl;


  /* NOTE
   * what if both left and right loci are under selection, potentially with different selection coefficients, even opposite signs
   * selection constrained to a particular epoch
   * epoch-specific mutation rates to test the hypothesis of change in mu along human lineage (wrt eg chimps)?
   * alias N's from adjacent epochs to reduce parameter space
   */

  if(argc == 1)
  {
    std::cout << "Please fill in a text file with the following options and execute from the command line:\nmomentspp params=[file]\n\n";

    std::cout << "label = # optional, default = 'moments++'\n";
    std::cout << "demes_file = # mandatory, relative path to file in Demes format that specifies the (starting) model\n";
    std::cout << "stats_file = # optional, relative path to file listing observed summary statistics from sampled populations, default = 'none'\n\n";

    std::cout << "optimizer = # optional, default = 'BFGS'\n";
    std::cout << "tolerance = # optional, default = 1e-6\n";
    std::cout << "compress_moments = # optional, default = 1 (TRUE)\n";
    std::cout << "num_threads = # optional, default = num_cores / 2\n";

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

    size_t start = demes.getPopsVec()[i].front()->getStartTime();
    size_t end = demes.getPopsVec()[i].front()->getEndTime();

    SumStatsLibrary sslib(demes.getPopsVec()[i], options.compressMoments());

    /* Epoch-specific operators (concern populations present in each epoch, hence parameters must follow suit)
     * must have epoch-specific recombination and mutation operators because they depend on pop indices (popss[i]),
     * even though inside Model we alias r and mu across epochs
     */

    std::vector<double> drift(0);
    drift.reserve(demes.getPopsVec()[i].size());

    // from (diploid) population sizes (N_j, not 2N_j) to drift parameters
    for(size_t j = 0; j < demes.getPopsVec()[i].size(); ++j)
      drift.emplace_back(1. / (2. * demes.getPopsVec()[i][j]->getSize()));

    std::shared_ptr<bpp::IntervalConstraint> ic = std::make_shared<bpp::IntervalConstraint>(0., 1e-2, true, true);
    std::shared_ptr<Drift> driftOp = std::make_shared<Drift>(drift, ic, sslib);
    std::shared_ptr<Recombination> recOp = std::make_shared<Recombination>(demes.getRec(i), ic, sslib);
    std::shared_ptr<Mutation> mutOp = std::make_shared<Mutation>(demes.getMu(i), ic, sslib);
    std::shared_ptr<Admixture> admixOp = nullptr;

    std::vector<std::shared_ptr<AbstractOperator>> operators(0);
    operators.reserve(4);

    if(demes.getNumPops(i) > 1)
    {
      std::shared_ptr<Migration> migOp = std::make_shared<Migration>(demes.getMig(i), ic, sslib);
      operators.emplace_back(migOp);

      if(!demes.getPulse(i).isZero(0))
      {
        admixOp = std::make_shared<Admixture>(demes.getPulse(i), sslib);
        //admixOp->printDeltaLDMat(id + "_admix.csv", sslib);
      }
    }

    operators.emplace_back(driftOp);
    operators.emplace_back(recOp);
    operators.emplace_back(mutOp);

    epochs.emplace_back(std::make_shared<Epoch>(id, sslib, start, end, operators, admixOp, demes.getPopsVec()[i]));
    epochs.back()->printAttributes(std::cout);
  }

  epochs[0]->pseudoSteadyState(); // only need to have steady state in the deepest epoch

  try
  {
    if(options.getDataFilePath() == "none")
    {
      std::cout << "\nNo stats_file provided, moments++ will output expectations for input parameters.\n";
      std::shared_ptr<Model> model = std::make_shared<Model>(options.getLabel(), epochs);

      //model->getParameters().printParameters(std::cout);
      //model->getIndependentParameters().printParameters(std::cout);
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
