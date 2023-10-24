/*
 * Author: Gustavo V. Barroso
 * Created: 20/10/2023
 * Last modified: 24/10/2023
 * Source code for twoLocusSim
 *
 */

////////////////////////////////////////////////////////////////////////////////////////////////////

#include <vector>
#include <map>
#include <limits>
#include <thread>

#include <stdio.h>
#include <sys/time.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Text/TextTools.h>

#include "TwoLocusPair.hpp"
#include "TwoLocusPop.hpp"


int main(int argc, char *argv[]) {

  std::cout << std::endl;
  std::cout << "******************************************************************" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "*                TwoLocusSim  version 0.0.1                      *" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "*            Two-site recursions                                 *" << std::endl;
  std::cout << "*            Unraveling history                                  *" << std::endl;
  std::cout << "*            Moment by moment                                    *" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "* Authors: G. V. Barroso                 Last Modif. 24/Oct/2023 *" << std::endl;
  std::cout << "*          A. P. Ragsdale                                        *" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "******************************************************************" << std::endl;

  std::cout << "\nCompiled on: " << __DATE__ << std::endl;
  std::cout << "Compiled at: " << __TIME__ << std::endl << std::endl;

  if(argc == 1)
  {
    std::cout << "To use TwoLocusSim, fill in a text file with the following options and execute from the command line:\ttwolocussim params=file_name\n\n";

    std::cout << "L = \n";
    std::cout << "N = \n";
    std::cout << "u = \n";
    std::cout << "r = \n";
    std::cout << "s = \n";

    std::cout << "\nFor more information, please email gvbarroso@gmail.com " << std::endl;
    return(0);
  }

  bpp::BppApplication twoLocusSim(argc, argv, "TwoLocusSim");
  twoLocusSim.startTimer();
  std::map<std::string, std::string> params = twoLocusSim.getParams();

  // TODO make N's and G's as vector parameters
  size_t G = bpp::ApplicationTools::getParameter<size_t>("G", params, 1000000, "", 0);
  size_t L = bpp::ApplicationTools::getParameter<size_t>("L", params, 1, "", 0);
  unsigned int N = bpp::ApplicationTools::getParameter<unsigned int>("N", params, 10000, "", 0);
  double u = bpp::ApplicationTools::getParameter<double>("u", params, 1e-6, "", 0);
  double r = bpp::ApplicationTools::getParameter<double>("r", params, 1e-7, "", 0);
  double s = bpp::ApplicationTools::getParameter<double>("s", params, -1e-4, "", 0);

  struct timeval tv;
  gettimeofday(&tv, 0);
  unsigned long seed = tv.tv_sec + tv.tv_usec;

  const gsl_rng_type * T;
  gsl_rng * gen;

  gsl_rng_env_setup();

  T = gsl_rng_default;
  gen = gsl_rng_alloc(T);

  gsl_rng_set(gen, seed);

  std::vector<TwoLocusPair> pairs(0);
  pairs.reserve(2 * N);

  TwoLocusPop pop(0, L, N, pairs);

  std::cout << "Burn-in...\n";

  for(size_t i = 0; i < 1e+6; ++i)
    pop.evolve_random(gen, i, u, r, s);

  std::cout << "done. \nPerforming simulation...\n";

  std::array<double, 5> gen_stats;
  double avg_Hl = 0.;
  double avg_Hr = 0.;
  double avg_pi2 = 0.;
  double avg_Dz = 0.;
  double avg_Dsqr = 0.;

  for(size_t i = 0; i < G; ++i)
  {
    pop.evolve_random(gen, i, u, r, s);
    gen_stats = pop.fetchAvgStats();

    avg_Hl += gen_stats[0];
    avg_Hr += gen_stats[1];
    avg_pi2 += gen_stats[2];
    avg_Dz += gen_stats[3];
    avg_Dsqr += gen_stats[4];
  }

  std::cout << "avg_Hl = " << avg_Hl / G << "\n" ;
  std::cout << "avg_Hr = " << avg_Hr / G << "\n" ;
  std::cout << "avg_pi2 = " << avg_pi2 / G << "\n" ;
  std::cout << "avg_Dz = " << avg_Dz / G << "\n" ;
  std::cout << "avg_Dsqr = " << avg_Dsqr / G << "\n" ;

  gsl_rng_free(gen);

  twoLocusSim.done();
  return 0;
}
