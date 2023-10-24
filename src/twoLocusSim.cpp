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

  // TODO make Population class, N's and G's as vector parameters
  size_t G = bpp::ApplicationTools::getParameter<size_t>("G", params, 1000000, "", 0);
  size_t L = bpp::ApplicationTools::getParameter<size_t>("L", params, 1000000, "", 0);
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
  pairs.reserve(1e+2 * N);

  for(size_t i = 0; i < G; ++i) // for each epoch, from past to present
  {
    size_t unlinkedMuts = gsl_ran_poisson(gen, L * N * u);

    for(size_t j = 0; j < unlinkedMuts; ++j)
    {
      unsigned int c_ab = N - 1;
      unsigned int c_Ab = 0;
      unsigned int c_aB = 0;
      unsigned int c_AB = 0;

      if(gsl_rng_uniform(gen) < 0.5)
        c_Ab = 1;

      else
        c_aB = 1;

      TwoLocusPair newPair(i, c_ab, c_Ab, c_aB, c_AB);
      pairs.emplace_back(newPair);
    }

    std::cout << "Gen. " << i << ": " << pairs.size() << " segregating pairs.\n" ;

    // summary statistics
    double sum_hl = 0.;
    double sum_hr = 0.;
    double sum_pi2 = 0.;
    double sum_dz = 0.;
    double sum_dsqr = 0.;

    for(auto it = std::begin(pairs); it != std::end(pairs);)
    {
      it->evolve_random(gen, u, r, s);
      //it->printAttributes(std::cout);

      sum_hl += it->fetchHl();
      sum_hr += it->fetchHr();

      if(it->bothPolymorphic())
      {
        sum_pi2 += it->fetchPi2();
        sum_dz += it->fetchDz();
        sum_dsqr += it->fetchDsqr();
      }

      if(it->monomorphic() && it->mutatedBoth())
        it = pairs.erase(it);

      else
        ++it;
    }

    /*std::cout << "avg Hl = " << sum_hl / pairs.size() << "\n";
    std::cout << "avg Hr = " << sum_hr / pairs.size() << "\n";
    std::cout << "avg pi2 = " << sum_pi2 / pairs.size() << "\n";
    std::cout << "avg Dz = " << sum_dz / pairs.size() << "\n";
    std::cout << "avg D^2 = " << sum_dsqr / pairs.size() << "\n";*/
  }

  gsl_rng_free(gen);

  twoLocusSim.done();
  return 0;
}
