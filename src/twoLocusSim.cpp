/*
 * Author: Gustavo V. Barroso
 * Created: 20/10/2023
 * Last modified: 17/11/2023
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

// Functions
double getD(const std::vector<double>& X)
{
  return X[0] * X[3] - X[1] * X[2];
}

double getP(const std::vector<double>& X)
{
  return X[0] + X[1];
}

double getQ(const std::vector<double>& X)
{
  return X[0] + X[2];
}

double getHl(double Xl)
{
  return Xl * (1. - Xl);
}

double getHr(double Xr)
{
  return Xr * (1. - Xr);
}

double getDsqr(const std::vector<double>& X)
{
  return (X[0] * X[3] - X[1] * X[2]) * (X[0] * X[3] - X[1] * X[2]);
}

double getDz(const std::vector<double>& X)
{
  return getD(X) * (1. - 2. * getP(X)) * (1. - 2. * getQ(X));
}

double getPi2(const std::vector<double>& X)
{
  double p = getP(X);
  double q = getQ(X);
  return p * (1. - p) * q * (1. - q);
}

void recombine(std::vector<double>& X, double r)
{
  double d = getD(X);

  X[0] -= r * d;
  X[3] -= r * d;
  X[1] += r * d;
  X[2] += r * d;
}

void select(double& Xl, double s)
{
  Xl = ((1 + s) * Xl) / (Xl * (1 + s) + 1 - Xl);
}

void drift(std::vector<double>& X, double& Xl, double& Xr, double Ne, const gsl_rng* gen)
{
  Xl = gsl_ran_binomial(gen, Xl, 2. * Ne) / (2. * Ne);
  Xr = gsl_ran_binomial(gen, Xr, 2. * Ne) / (2. * Ne);

  unsigned int n = std::round(std::accumulate(std::begin(X), std::end(X), 0.));
  unsigned int next[4];
  double probs[4] = X;

  gsl_ran_multinomial(gen, 4, n, probs, next);

  X[0] = next[0];
  X[1] = next[1];
  X[2] = next[2];
  X[3] = next[3];
}

void mutate()
{
}

void evolve(double Xl, double Xr, std::vector<double>& X, double Ne, double u, double r, double L, double s, const gsl_rng* gen)
{
  recombine(X, r);
  select(Xl, s);
  drift(X, Xl, Xr, Ne, gen);
}

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
  std::cout << "* Authors: G. V. Barroso                 Last Modif. 09/Nov/2023 *" << std::endl;
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
  size_t B = bpp::ApplicationTools::getParameter<size_t>("B", params, 1, "", 0);
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
  pairs.reserve(2 * L * N * u);

  TwoLocusPop root(0, L, N, pairs);

  std::cout << "\nBurn-in (" << B << " generations)...";
  std::cout.flush();

  for(size_t i = 0; i < B; ++i)
    root.evolve_det(gen, u, r, s);

  std::cout << "done.\n";

  TwoLocusPop p1 = root;
  //TwoLocusPop p2 = root;

  //p1.setPopSize(N * 2);
  //p2.setPopSize(N * 10);

  p1.setId(1);
  //p2.setId(2);

  std::cout << "\nEvolving population(s) (" << G << " generations)...\n";

  std::array<double, 5> gen_stats_1;
  double avg_Hl_1 = 0.;
  double avg_Hr_1 = 0.;
  double avg_pi2_1 = 0.;
  double avg_Dz_1 = 0.;
  double avg_Dsqr_1 = 0.;

  /*std::array<double, 5> gen_stats_2;
  double avg_Hl_2 = 0.;
  double avg_Hr_2 = 0.;
  double avg_pi2_2 = 0.;
  double avg_Dz_2 = 0.;
  double avg_Dsqr_2 = 0.;*/

  for(size_t i = 0; i < G; ++i)
  {
    p1.evolve_det(gen, u, r, s);
    //p2.evolve_det(gen, i, u, r, s);

    gen_stats_1 = p1.fetchAvgStats();
    //gen_stats_2 = p2.fetchAvgStats();

    avg_Hl_1 += gen_stats_1[0];
    avg_Hr_1 += gen_stats_1[1];

    // NOTE divide two-locus stats by Nu to bring them to the same scale as H's.
    // This is required because of the nature of the infinite sites model.
    // In this simulator, each new mutation creates a two-locus system,
    // but at this time only one of the loci had the opportunity to mutate.
    // Meanwhile, two-locus stats are still computed and averaged over generations
    avg_pi2_1 += gen_stats_1[2] / (N * u);
    avg_Dz_1 += gen_stats_1[3] / (N * u);
    avg_Dsqr_1 += gen_stats_1[4] / (N * u);

    /*avg_Hl_2 += gen_stats_2[0];
    avg_Hr_2 += gen_stats_2[1];
    avg_pi2_2 += gen_stats_2[2];
    avg_Dz_2 += gen_stats_2[3];
    avg_Dsqr_2 += gen_stats_2[4];*/
  }

  std::cout << "avg_Hl_1 = " << avg_Hl_1 / G << "\n" ;
  std::cout << "avg_Hr_1 = " << avg_Hr_1 / G << "\n" ;
  std::cout << "avg_pi2_1 = " << avg_pi2_1 / G << "\n" ;
  std::cout << "avg_Dz_1 = " << avg_Dz_1 / G << "\n" ;
  std::cout << "avg_Dsqr_1 = " << avg_Dsqr_1 / G << "\n\n" ;

  /*std::cout << "avg_Hl_2 = " << avg_Hl_2 / G << "\n" ;
  std::cout << "avg_Hr_2 = " << avg_Hr_2 / G << "\n" ;
  std::cout << "avg_pi2_2 = " << avg_pi2_2 / G << "\n" ;
  std::cout << "avg_Dz_2 = " << avg_Dz_2 / G << "\n" ;
  std::cout << "avg_Dsqr_2 = " << avg_Dsqr_2 / G << "\n" ;*/

  gsl_rng_free(gen);

  twoLocusSim.done();
  return 0;
}
