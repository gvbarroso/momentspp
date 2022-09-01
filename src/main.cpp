/*
 * Author: Gustavo V. Barroso
 * Created: 29/08/2022
 * Last modified: 31/08/2022
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
  std::cout << "*                moments++  version 0.0.1                        *" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "*            Recombination                                       *" << std::endl;
  std::cout << "*            A mosaic in the genome                              *" << std::endl;
  std::cout << "*            Endless ancestors                                   *" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "* Authors: G. Barroso                    Last Modif. 01/Sep/2022 *" << std::endl;
  std::cout << "*          A. Ragsdale                                           *" << std::endl;
  std::cout << "******************************************************************" << std::endl;
  std::cout << std::endl;

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

    SumStatsLibrary ssl;
    ssl.init(data);

    std::cout << "done." << std::endl;

    OptimizationWrapper optimizer(options);

    if(options.resume()) // if resume optimization, read "backup_params.txt" from current dir
      optimizer.resumeOptim(ssl);

    else
      optimizer.optimize(ssl); // else optimize from scratch, backing up to "backup_params.txt"
  }

  catch(std::exception& e)
  {
    std::cout << "moments++ terminated because of an error." << std::endl;
    std::cout << e.what() << std::endl;
    return 1;
  }

  momentspp.done();
  return 0;
}
