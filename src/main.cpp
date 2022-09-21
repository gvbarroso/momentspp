/*
 * Author: Gustavo V. Barroso
 * Created: 29/08/2022
 * Last modified: 21/09/2022
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
  std::cout << "*            Two-locus histories                                 *" << std::endl;
  std::cout << "*            Revealed recursively                                *" << std::endl;
  std::cout << "*            Moment by moment                                    *" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "* Authors: G. Barroso                    Last Modif. 21/Sep/2022 *" << std::endl;
  std::cout << "*          A. Ragsdale                                           *" << std::endl;
  std::cout << "******************************************************************" << std::endl;
  std::cout << std::endl;

  /* TODO
  * 1. make it such (and ensure) that in every Epoch populations aren indexed from 0 to numPops - 1 (rendering popIndices_ obsolete)
  * This way we can avoid the binary search inside SumStatsLibrary setMomentValue methods (using pop id's directly instead)\
  * (also std::map<size_t, std::shared_ptr<Population>> pops_ in Epoch)
  */

  if(argc == 1)
  {
    std::cout << "To use moments++, please fill in the params file and simply call it from the command line: moments++ params=[params_file].bpp" << std::endl;
    std::cout << "For more information, please email gvbarroso@gmail.com " << std::endl;
    return(0);
  }

  BppApplication momentspp(argc, argv, "moments++"); // params file holding users choices (see OptionsContainer class)
  momentspp.startTimer();
  std::map<std::string, std::string> params = momentspp.getParams();

  OptionsContainer options(params);

  try
  {
    std::cout << "Processing input data..."; std::cout.flus();

    PolymorphismData data; // input data, format = ?
    data.parse(options.getDataPath());
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
