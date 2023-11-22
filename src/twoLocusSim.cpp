/*
 * Author: Gustavo V. Barroso
 * Created: 20/10/2023
 * Last modified: 22/11/2023
 * Source code for twoLocusSim
 *
 */

////////////////////////////////////////////////////////////////////////////////////////////////////

#include <vector>
#include <limits>
#include <thread>
#include <array>
#include <numeric>

#include <stdio.h>
#include <sys/time.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>

#include "TwoLocusPop.hpp"

int main(int argc, char *argv[]) {

  std::cout << std::endl;
  std::cout << "******************************************************************" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "*                TwoLocusSim  version 0.0.1                      *" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "*            Two-site simulations                                *" << std::endl;
  std::cout << "*            A companion for moments++                           *" << std::endl;
  std::cout << "*            This is not a haiku                                 *" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "* Authors: G. V. Barroso                 Last Modif. 22/Nov/2023 *" << std::endl;
  std::cout << "*          A. P. Ragsdale                                        *" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "******************************************************************" << std::endl;

  std::cout << "\nCompiled on: " << __DATE__ << std::endl;
  std::cout << "Compiled at: " << __TIME__ << std::endl << std::endl;

  if(argc == 1)
  {
    std::cout << "To use TwoLocusSim, fill in a text file with the following options and execute from the command line:\ttwolocussim params=file_name\n\n";

    std::cout << "L = \n";
    std::cout << "Ne = \n";
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
  unsigned int Ne = bpp::ApplicationTools::getParameter<unsigned int>("Ne", params, 10000, "", 0);
  double u = bpp::ApplicationTools::getParameter<double>("u", params, 1e-6, "", 0);
  double r = bpp::ApplicationTools::getParameter<double>("r", params, 1e-7, "", 0);
  double s = bpp::ApplicationTools::getParameter<double>("s", params, 0., "", 0);

  unsigned int B = 40 * Ne;

  struct timeval tv;
  gettimeofday(&tv, 0);
  unsigned long seed = tv.tv_sec + tv.tv_usec;

  const gsl_rng_type * T;
  gsl_rng * gen;

  gsl_rng_env_setup();

  T = gsl_rng_default;
  gen = gsl_rng_alloc(T);

  gsl_rng_set(gen, seed);

  TwoLocusPop root(0, L, Ne);

  std::cout << "\nBurn-in (" << B << " generations)...";
  std::cout.flush();

  for(size_t i = 0; i < B; i++)
  {
    root.evolve(u, r, s, gen);
    root.cleanup();
  }

  TwoLocusPop p1 = root;
  //TwoLocusPop p2 = root;

  std::cout << "done.\n";
  std::cout << "\nEvolving population(s) (" << G << " generations)...\n";

  // TODO multi-thread from root in 2-pops case

  double sum_Hl = 0.;
  double sum_Hr = 0.;
  double sum_pi2 = 0.;
  double sum_Dz = 0.;
  double sum_Dsqr = 0.;

  for(size_t i = 0; i < G; ++i)
  {
    p1.evolve(u, r, s, gen);
    p1.cleanup();

    sum_Hl += p1.getSumHl();
    sum_Hr += p1.getSumHr();
    sum_Dz += p1.getSumDz();
    sum_Dsqr += p1.getSumDsqr();
    sum_pi2 += p1.getSumPi2();

    if(i % (10 * Ne) == 1)
    {
      std::cout << "Generation " << i << "\n";
      std::cout << "\tavg_Hl = " << sum_Hl / L / i << "\n";
      std::cout << "\tavg_Hr = " << sum_Hr / L / i << "\n";
      std::cout << "\tavg_Dz = " << sum_Dz / L / L / i << "\n";
      std::cout << "\tavg_Dsqr = " << sum_Dsqr / L / L / i << "\n";
      std::cout << "\tavg_pi2 = " << sum_pi2 / L / L / i << "\n\n";
    }
  }

  std::cout << "avg_Hl = " << sum_Hl / L / G << "\n";
  std::cout << "avg_Hr = " << sum_Hr / L / G << "\n";
  std::cout << "avg_Dz = " << sum_Dz / L / L / G << "\n";
  std::cout << "avg_Dsqr = " << sum_Dsqr / L / L / G << "\n";
  std::cout << "avg_pi2 = " << sum_pi2 / L / L / G << "\n\n";

  gsl_rng_free(gen);

  twoLocusSim.done();
  return 0;
}
