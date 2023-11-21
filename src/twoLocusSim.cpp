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
#include <array>

#include <stdio.h>
#include <sys/time.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Text/TextTools.h>

// Functions
std::vector<double> getDs(const std::vector<std::array<double, 4>>& X)
{
  std::vector<double> ret(0);
  ret.reserve(X.size());

  for(size_t i = 0; i < X.size(); ++i)
    ret.emplace_back(X[i][0] * X[i][3] - X[i][1] * X[i][2]);

  return ret;
}

std::vector<double> getPs(const std::vector<std::array<double, 4>>& X)
{
  std::vector<double> ret(0);
  ret.reserve(X.size());

  for(size_t i = 0; i < X.size(); ++i)
    ret.emplace_back(X[i][0] + X[i][1]);

  return ret;
}

std::vector<double> getQs(const std::vector<std::array<double, 4>>& X)
{
  std::vector<double> ret(0);
  ret.reserve(X.size());

  for(size_t i = 0; i < X.size(); ++i)
    ret.emplace_back(X[i][0] + X[i][2]);

  return ret;
}

std::vector<double> getHls(const std::vector<double>& Xl)
{
  std::vector<double> ret(0);
  ret.reserve(Xl.size());

  for(size_t i = 0; i < Xl.size(); ++i)
    ret.emplace_back(Xl[i] * (1. - Xl[i]));

  return ret;
}

std::vector<double> getHrs(const std::vector<double>& Xr)
{
  std::vector<double> ret(0);
  ret.reserve(Xr.size());

  for(size_t i = 0; i < Xr.size(); ++i)
    ret.emplace_back(Xr[i] * (1. - Xr[i]));

  return ret;
}

std::vector<double> getDsqrs(const std::vector<std::array<double, 4>>& X)
{
  std::vector<double> ret(0);
  ret.reserve(X.size());

  for(size_t i = 0; i < X.size(); ++i)
    ret.emplace_back((X[i][0] * X[i][3] - X[i][1] * X[i][2]) * (X[i][0] * X[i][3] - X[i][1] * X[i][2]));

  return ret;
}

std::vector<double> getDzs(const std::vector<std::array<double, 4>>& X)
{
  std::vector<double> ret(0);
  ret.reserve(X.size());

  std::vector<double> ps = getPs(X);
  std::vector<double> qs = getQs(X);
  std::vector<double> ds = getDs(X);

  for(size_t i = 0; i < X.size(); ++i)
    ret.emplace_back(ds[i] * (1. - 2. * ps[i]) * (1. - 2. * qs[i]));

  return ret;
}

std::vector<double> getPi2s(const std::vector<std::array<double, 4>>& X)
{
  std::vector<double> ret(0);
  ret.reserve(X.size());

  std::vector<double> ps = getPs(X);
  std::vector<double> qs = getQs(X);

  for(size_t i = 0; i < X.size(); ++i)
    ret.emplace_back(ps[i] * (1. - ps[i]) * qs[i] * (1. - qs[i]));

  return ret;
}

void recombine(std::vector<std::array<double, 4>>& X, double r)
{
  std::vector<double> ds = getDs(X);

  for(size_t i = 0; i < X.size(); ++i)
  {
    X[i][0] -= r * ds[i];
    X[i][3] -= r * ds[i];
    X[i][1] += r * ds[i];
    X[i][2] += r * ds[i];
  }
}

void select(std::vector<double>& Xl, double s)
{
  for(size_t i = 0; i < Xl.size(); ++i)
    Xl[i] = ((1 + s) * Xl[i]) / (Xl[i] * (1 + s) + 1 - Xl[i]);
}

void drift(std::vector<std::array<double, 4>>& X, std::vector<double>& Xl, std::vector<double>& Xr, double Ne, const gsl_rng* gen)
{
  for(size_t i = 0; i < Xl.size(); ++i)
    Xl[i] = gsl_ran_binomial(gen, Xl[i], 2. * Ne) / (2. * Ne);

  for(size_t i = 0; i < Xr.size(); ++i)
    Xr[i] = gsl_ran_binomial(gen, Xr[i], 2. * Ne) / (2. * Ne);

  unsigned int n = std::round(2. * Ne);

  for(size_t i = 0; i < X.size(); ++i)
  {
    unsigned int next[4];
    double probs[4] = X[i]; // NOTE can be converted?

    gsl_ran_multinomial(gen, 4, n, probs, next);

    X[i][0] = next[0];
    X[i][1] = next[1];
    X[i][2] = next[2];
    X[i][3] = next[3];
  }
}

void mutate(std::vector<std::array<double, 4>& X, std::vector<double>& Xl, std::vector<double>& Xr, double Ne, double u, double L, const gsl_rng* gen)
{
  // new single mutations at left locus
  unsigned int num_left_mut = gsl_ran_poisson(2 * Ne * u * L);
  for(size_t i = 0; i < num_left_mut; ++i)
    Xl.emplace_back(1. / (2. * Ne));

  // new single mutations at right locus
  unsigned int num_right_mut = gsl_ran_poisson(2 * Ne * u * L);
  for(size_t i = 0; i < num_right_mut; ++i)
    Xr.emplace_back(1. / (2. * Ne));


  for(size_t i = 0; i < L; ++i)
  {
    // mutate against segregating left loci in Xl
    for(size_t j = 0; j < Xl.size(); ++j)
    {
      if(gsl_rng_uniform(gen) < 2 * Ne * u)
      {
        if(gsl_rng_uniform(gen) < Xl[j]) // falls on bA background
          X.push_back(1. / 2 / Ne, Xl[j] - 1 / 2 / Ne, 0, 1 - Xl[j]]])
                    )
        else // falls on ab background
          X = np.concatenate(
                        (X, [[0, Xl[j], 1 / 2 / Ne, 1 - Xl[j] - 1 / 2 / Ne]])
                    )
      }
    }

    // mutate against segregating right loci in Xr
    for(size_t j = 0; j < Xr.size(); ++j)
    {
      if(gsl_rng_uniform(gen) < 2 * Ne * u)
      {
        if(gsl_rng_uniform(gen) < Xr[j]) // falls on bA background
        // falls on aB background
            X = , [[1 / 2 / Ne, 0, Xr[j] - 1 / 2 / Ne, 1 - Xr[j]]])
                    )
                else:
                    // falls on ab background
                    X = np.concatenate(
                        (X, [[0, 1 / 2 / Ne, Xr[j], 1 - Xr[j] - 1 / 2 / Ne]])
                    )

      }
    }
  }
}

void evolve(std::vector<double>& Xl, std::vector<double>& Xr, std::vector<double>& X, double Ne, double u, double r, double L, double s, const gsl_rng* gen)
{
  recombine(X, r);
  select(Xl, s);
  drift(X, Xl, Xr, Ne, gen);
  mutate(X, Xl, Xr, Ne, u, L, gen);
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

  std::array<double, 5> gen_stats_2;
  double avg_Hl_2 = 0.;
  double avg_Hr_2 = 0.;
  double avg_pi2_2 = 0.;
  double avg_Dz_2 = 0.;
  double avg_Dsqr_2 = 0.;

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

    avg_Hl_2 += gen_stats_2[0];
    avg_Hr_2 += gen_stats_2[1];
    avg_pi2_2 += gen_stats_2[2];
    avg_Dz_2 += gen_stats_2[3];
    avg_Dsqr_2 += gen_stats_2[4];
  }

  std::cout << "avg_Hl_1 = " << avg_Hl_1 / G << "\n" ;
  std::cout << "avg_Hr_1 = " << avg_Hr_1 / G << "\n" ;
  std::cout << "avg_pi2_1 = " << avg_pi2_1 / G << "\n" ;
  std::cout << "avg_Dz_1 = " << avg_Dz_1 / G << "\n" ;
  std::cout << "avg_Dsqr_1 = " << avg_Dsqr_1 / G << "\n\n" ;

  std::cout << "avg_Hl_2 = " << avg_Hl_2 / G << "\n" ;
  std::cout << "avg_Hr_2 = " << avg_Hr_2 / G << "\n" ;
  std::cout << "avg_pi2_2 = " << avg_pi2_2 / G << "\n" ;
  std::cout << "avg_Dz_2 = " << avg_Dz_2 / G << "\n" ;
  std::cout << "avg_Dsqr_2 = " << avg_Dsqr_2 / G << "\n" ;

  gsl_rng_free(gen);

  twoLocusSim.done();
  return 0;
}
