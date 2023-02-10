/*
 * Author: Gustavo V. Barroso
 * Created: 29/08/2022
 * Last modified: 01/02/2023
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
  std::cout << "* Authors: G. Barroso                    Last Modif. 09/Feb/2023 *" << std::endl;
  std::cout << "*          A. Ragsdale                                           *" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "******************************************************************" << std::endl;

  std::cout << "\nCompiled on: " << __DATE__ << std::endl;
  std::cout << "Compiled at: " << __TIME__ << std::endl << std::endl;


  /* NOTE s
  * Write Admixture operator as an AbstractOperator that has exponent 1!
  * Selection operator -> upon rejection, sample a replacement individual from the whole population with probaility proportional to the pop. fitness vector (which maybe can be obtained from allele frequencies?)
  *
  * selection constrained to a particular epoch
  *
  * define start and end of epochs as quantiles of the exp dist?
  *
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
      bool hasSelection = 0; //j == 0;

      //Population(const std::string& name, const std::string& description, size_t id, size_t startTime, size_t endTime, size_t startSize, size_t endSize, bool hasSelection):
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

  try
  {
    if(options.getDataFilePath() == "none") // default
    {
      // just compute moments for models specified in demes files
    }

    else // there is a data file with observed summary statistics
    {
      Data data(options.getDataFilePath(), popMaps.back()); // input summary statistics (observed), format = ?

      // the optimizer builds the main objects, assembles the Models and optimizes them
      OptimizationWrapper optimizer(options);
      optimizer.optimize(data, demes);
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
