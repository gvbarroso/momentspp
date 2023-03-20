/*
 * Author: Gustavo V. Barroso
 * Created: 29/08/2022
 * Last modified: 20/03/2023
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
//#include "Admixture.hpp"
#include "OptimizationWrapper.hpp"
#include "OptionsContainer.hpp"
#include "Model.hpp"
#include "Data.hpp"
#include "Demes.hpp"

int main(int argc, char *argv[]) {

  std::cout << std::endl;
  std::cout << "******************************************************************" << std::endl;
  std::cout << "*                moments++  version 0.0.1                        *" << std::endl;
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
  std::cout << "* Authors: G. Barroso                    Last Modif. 20/Mar/2023 *" << std::endl;
  std::cout << "*          A. Ragsdale                                           *" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "******************************************************************" << std::endl;

  std::cout << "\nCompiled on: " << __DATE__ << std::endl;
  std::cout << "Compiled at: " << __TIME__ << std::endl << std::endl;


  /* NOTE
   * IDEAS:
   * Selection operator -> upon rejection, sample a replacement individual from the whole population with probaility proportional to the pop. fitness vector (which maybe can be obtained from allele frequencies?)
   * selection constrained to a particular epoch
   * define start and end of epochs as quantiles of the exp dist?
   *
   * TODO
   * Write Admixture operator as an AbstractOperator that has exponent 1
   * > 2 pops
   */

  if(argc == 1)
  {
    std::cout << "To use moments++, please fill in the params file and simply call it from the command line: momentspp params=[params_file].bpp" << std::endl;
    std::cout << "For more information, please email gvbarroso@gmail.com " << std::endl;
    return(0);
  }

  bpp::BppApplication momentspp(argc, argv, "moments++"); // params file holding users choices (see OptionsContainer class)
  momentspp.startTimer();
  std::map<std::string, std::string> params = momentspp.getParams();

  OptionsContainer options(params);

  // 1. parse options.getPopsFilePath()
  // 2. create populations
  // 3. link populations from different epochs (see Model::linkMoments())
  size_t numEpochs = options.getNumEpochs();
  std::vector<std::map<size_t, std::shared_ptr<Population>>> popMaps(0); // pop_id -> Population*, one per epoch
  popMaps.reserve(numEpochs);

  for(size_t i = 0; i < numEpochs; ++i)
  {
    std::map<size_t, std::shared_ptr<Population>> map;

    for(size_t j = 0; j < options.getNumPops(); ++j) // simplification: for now, every epoch has same number of populations; use Demes class to change that
    {
      bool hasSelection = 0; //j % 2 == 0;

      std::shared_ptr<Population> pop = std::make_shared<Population>("pop_" + bpp::TextTools::toString(j), "test population", j, 500000, 0, 10000, 10000, hasSelection);
      map.try_emplace(j, pop);
    }

    if(i > 0)
    {
      for(size_t j = 0; j < options.getNumPops(); ++j)
      {
        map.at(j)->setLeftParent(popMaps.back().at(j));
        map.at(j)->setRightParent(popMaps.back().at(j));
      }
    }

    popMaps.emplace_back(map);
  }

  Demes demes(popMaps, options.getDemesFilePath());

  std::vector<std::shared_ptr<Epoch>> epochs(0);
  epochs.reserve(numEpochs);

  std::shared_ptr<bpp::IntervalConstraint> ic = std::make_shared<bpp::IntervalConstraint>(0., 1., true, true); // NOTE

  for(size_t i = 0; i < numEpochs; ++i) // for each epoch, from past to present
  {
    std::string id = "e_" + bpp::TextTools::toString(i);

    size_t start = i * (options.getTotalNumberOfGenerations() / numEpochs);
    size_t end = (i + 1) * (options.getTotalNumberOfGenerations() / numEpochs);

    SumStatsLibrary sslib(options.getOrder(), demes.getPopMaps()[i]);

    /* Epoch-specific operators (concern populations present in each epoch, hence parameters must follow suit)
     * must have epoch-specific recombination and mutation operators because they depend on pop indices (popMaps[i]),
     * even though we prob. want single r and mu params in Model --> alias r and mu across epochs?
     */

    std::vector<double> drift = options.getInitPopSizes(); // should get this from Demes instead
    for(size_t j = 0; j < drift.size(); ++j) // from (diploid) population sizes (N_j, not 2N_j) to (diploid) drift parameters
      drift[j] = 1. / (2. * drift[j]);

    std::shared_ptr<Drift> driftOp = std::make_shared<Drift>(drift, ic, sslib);
    std::shared_ptr<Recombination> recOp = std::make_shared<Recombination>(options.getInitR(), ic, sslib);
    std::shared_ptr<Mutation> mutOp = std::make_shared<Mutation>(options.getInitMu(), ic, sslib);

    std::vector<std::shared_ptr<AbstractOperator>> operators(0);
    operators.reserve(4);

    if(demes.getNumPops() > 1)
    {
      std::shared_ptr<Migration> migOp = std::make_shared<Migration>(options.getInitMig(), ic, sslib);
      operators.emplace_back(migOp);
    }

    operators.emplace_back(driftOp);
    operators.emplace_back(recOp);
    operators.emplace_back(mutOp);

    epochs.emplace_back(std::make_shared<Epoch>(id, sslib, start, end, operators, demes.getPopMaps()[i]));

    #ifdef VERBOSE
    std::ofstream recOut;
    recOut.open(epochs.back()->getName() + "_recursions.txt");
    epochs.back()->printRecursions(recOut);
    recOut.close();
    #endif
  }

  try
  {
    if(options.getDataFilePath() == "none") // no data file (default), we just compute expectations for given input parameters
    {
      std::shared_ptr<Model> model = std::make_shared<Model>(options.getLabel(), epochs);
      model->computeExpectedSumStats();
      std::ofstream fout(model->getName() + "_final_unsorted.txt");
      model->printAliasedMoments(fout);
      fout.close();
    }

    else // there is a data file with observed summary statistics, thus we optimize
    {
      std::shared_ptr<Data> data = std::make_shared<Data>(options.getDataFilePath(), popMaps.back());
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
