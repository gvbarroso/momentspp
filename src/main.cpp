/*
 * Author: Gustavo V. Barroso
 * Created: 29/08/2022
 * Last modified: 23/09/2022
 * Source code for moments++
 *
 */

////////////////////////////////////////////////////////////////////////////////////////////////////

#include "SumStatsLibrary.hpp"
#include "Mutation.hpp"
#include "Recombination.hpp"
#include "Drift.hpp"
#include "Migration.hpp"
//#include "Selection.h"
//#include "Admixture.h"
#include "OptimizationWrapper.hpp"
#include "OptionsContainer.hpp"
#include "Model.hpp"
#include "PolymorphismData.hpp"

int main(int argc, char *argv[]) {

  std::cout << std::endl;
  std::cout << "******************************************************************" << std::endl;
  std::cout << "*                moments++  version 0.0.1                        *" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "*            Two-site recursions                                 *" << std::endl;
  std::cout << "*            Unraveling history                                  *" << std::endl;
  std::cout << "*            Moment by moment                                    *" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "* Authors: G. Barroso                    Last Modif. 27/Sep/2022 *" << std::endl;
  std::cout << "*          A. Ragsdale                                           *" << std::endl;
  std::cout << "******************************************************************" << std::endl;
  std::cout << std::endl;

  /* TODO
  * Create Derived classes from Moment (DD, Dz, H, Pi2)?
  * Write Admixture operator as an AbstractOperator that has exponent 1!
  * Selection operator -> upon rejection, sample a replacement individual from the whole population with probaility proportional to the pop. fitness vector (which maybe can be obtained from allele frequencies?)
  * Linked selection -> another angle at iSMC's results?
  */

  if(argc == 1)
  {
    std::cout << "To use moments++, please fill in the params file and simply call it from the command line: moments++ params=[params_file].bpp" << std::endl;
    std::cout << "For more information, please email gvbarroso@gmail.com " << std::endl;
    return(0);
  }

  bpp::BppApplication momentspp(argc, argv, "moments++"); // params file holding users choices (see OptionsContainer class)
  momentspp.startTimer();
  std::map<std::string, std::string> params = momentspp.getParams();

  OptionsContainer options(params);

  std::vector<std::map<size_t, std::shared_ptr<Population>>> popMaps(0); // one per epoch
  // 1. parse options.getPopsFilePath()
  // 2. create populations
  // 3. link populations (see Model::linkMoments())

  try
  {
    std::cout << "Processing input data..."; std::cout.flush();

    PolymorphismData data; // input data, format = ?
    data.parse(options.getDataFilePath());
    data.computeSumStats();

    std::cout << "done." << std::endl;

    OptimizationWrapper optimizer(options);
    optimizer.optimize(data); // else optimize from scratch, backing up to "backup_params.txt"
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
