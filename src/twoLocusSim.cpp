/*
 * Author: Gustavo V. Barroso
 * Created: 20/10/2023
 * Last modified: 21/11/2023
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
#include <Bpp/Text/TextTools.h>

// Functions
double getP(const std::array<double, 4>& haps)
{
  return haps[0] + haps[1];
}

double getQ(const std::array<double, 4>& haps)
{
  return haps[0] + haps[2];
}

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

double getSumHl(const std::vector<double>& Xl)
{
  std::vector<double> vals = getHls(Xl);
  std::sort(std::begin(vals), std::end(vals));
  return std::accumulate(std::begin(vals), std::end(vals), 0.);
}

double getSumHr(const std::vector<double>& Xr)
{
  std::vector<double> vals = getHrs(Xr);
  std::sort(std::begin(vals), std::end(vals));
  return std::accumulate(std::begin(vals), std::end(vals), 0.);
}

double getSumDsqr(const std::vector<std::array<double, 4>>& X)
{
  std::vector<double> vals = getDsqrs(X);
  std::sort(std::begin(vals), std::end(vals));
  return std::accumulate(std::begin(vals), std::end(vals), 0.);
}

double getSumDz(const std::vector<std::array<double, 4>>& X)
{
  std::vector<double> vals = getDzs(X);
  std::sort(std::begin(vals), std::end(vals));
  return std::accumulate(std::begin(vals), std::end(vals), 0.);
}

double getSumPi2(const std::vector<std::array<double, 4>>& X)
{
  std::vector<double> vals = getPi2s(X);
  std::sort(std::begin(vals), std::end(vals));
  return std::accumulate(std::begin(vals), std::end(vals), 0.);
}

void recombine(std::vector<std::array<double, 4>>& X, double r)
{
  std::vector<double> ds = getDs(X);

  for(size_t i = 0; i < X.size(); ++i)
  {
    X[i][0] -= r * ds[i];
    X[i][1] += r * ds[i];
    X[i][2] += r * ds[i];
    X[i][3] -= r * ds[i];
  }
}

void select(std::vector<double>& Xl, double s)
{
  for(size_t i = 0; i < Xl.size(); ++i)
    Xl[i] = ((1 + s) * Xl[i]) / (Xl[i] * (1 + s) + 1 - Xl[i]);
}

void drift(std::vector<std::array<double, 4>>& X, std::vector<double>& Xl, std::vector<double>& Xr, unsigned int Ne, const gsl_rng* gen)
{
  for(size_t i = 0; i < Xl.size(); ++i)
    Xl[i] = gsl_ran_binomial(gen, Xl[i], 2. * Ne) / (2. * Ne);

  for(size_t i = 0; i < Xr.size(); ++i)
    Xr[i] = gsl_ran_binomial(gen, Xr[i], 2. * Ne) / (2. * Ne);

  unsigned int n = std::round(2. * Ne);

  for(size_t i = 0; i < X.size(); ++i)
  {
    unsigned int next[4];
    double probs[4] = { X[i][0], X[i][1], X[i][2], X[i][3] };

    gsl_ran_multinomial(gen, 4, n, probs, next);

    X[i][0] = static_cast<double>(next[0]) / (2. * Ne);
    X[i][1] = static_cast<double>(next[1]) / (2. * Ne);
    X[i][2] = static_cast<double>(next[2]) / (2. * Ne);
    X[i][3] = static_cast<double>(next[3]) / (2. * Ne);
  }
}

void mutate(std::vector<std::array<double, 4>>& X, std::vector<double>& Xl, std::vector<double>& Xr, unsigned int Ne, size_t L, double u, const gsl_rng* gen)
{
  // new single mutations at left locus
  unsigned int num_left_mut = gsl_ran_poisson(gen, 2. * Ne * u * L);
  for(size_t i = 0; i < num_left_mut; ++i)
    Xl.emplace_back(1. / (2. * Ne));

  // new single mutations at right locus
  unsigned int num_right_mut = gsl_ran_poisson(gen, 2. * Ne * u * L);
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
        {
          std::array<double, 4> vals = {1. / (2. * Ne), Xl[j] - 1. / (2. * Ne), 0., 1. - Xl[j]};
          X.emplace_back(vals);
        }

        else // falls on ab background
        {
          std::array<double, 4> vals = {0., Xl[j], 1. / (2. * Ne), 1. - Xl[j] - 1. / (2. * Ne)};
          X.emplace_back(vals);
        }
      }
    }

    // mutate against segregating right loci in Xr
    for(size_t j = 0; j < Xr.size(); ++j)
    {
      if(gsl_rng_uniform(gen) < 2 * Ne * u)
      {
        if(gsl_rng_uniform(gen) < Xr[j]) // falls on bA background
        {
          std::array<double, 4> vals = {1. / (2. * Ne), Xr[j] - 1. / (2. * Ne), 0., 1. - Xr[j]};
          X.emplace_back(vals);
        }

        else // falls on ab background
        {
          std::array<double, 4> vals = {0., 1. / (2. * Ne), Xr[j], 1. - Xr[j] - 1. / (2. * Ne)};
          X.emplace_back(vals);
        }
      }
    }
  }
}

void evolve(std::vector<std::array<double, 4>>& X, std::vector<double>& Xl, std::vector<double>& Xr, unsigned int Ne, size_t L, double u, double r, double s, const gsl_rng* gen)
{
  recombine(X, r);
  select(Xl, s);
  drift(X, Xl, Xr, Ne, gen);
  mutate(X, Xl, Xr, Ne, L, u, gen);
}

void cleanup(std::vector<std::array<double, 4>>& X, std::vector<double>& Xl, std::vector<double>& Xr)
{
  // remove fixed or lost in Xl and Xr
  for(auto it = std::begin(Xl); it != std::end(Xl);)
  {
    if((*it) == 0. || (*it) == 1.)
      it = Xl.erase(it);

    else
      ++it;
  }

  for(auto it = std::begin(Xr); it != std::end(Xr);)
  {
    if((*it) == 0. || (*it) == 1.)
      it = Xr.erase(it);

    else
      ++it;
  }

  // remove any with p or q at 0 or 1
  for(auto it = std::begin(X); it != std::end(X);)
  {
    double p = getP(*it);
    double q = getQ(*it);

    if(p == 0. || p == 1. || q == 0. || q == 1.)
      it = X.erase(it);

    else
      ++it;
  }
}

void printXl(const std::vector<double>& Xl)
{
  std::cout << "Xl:\n\t";

  for(auto& v : Xl)
    std::cout << v << ",";

  std::cout << "\n";
}

void printXr(const std::vector<double>& Xr)
{
  std::cout << "Xr:\n\t";

  for(auto& v : Xr)
    std::cout << v << ",";

  std::cout << "\n";
}

void printX(const std::vector<std::array<double, 4>>& X)
{
  for(size_t i = 0; i < X.size(); ++i)
  {
    std::cout << "\tfAB = " << X[i][0] << ", ";
    std::cout << "fAb = " << X[i][1] << ", ";
    std::cout << "faB = " << X[i][2] << ", ";
    std::cout << "fab = " << X[i][3] << "\n";
  }

  std::cout << "\n";
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

  std::vector<double> Xl(0);
  std::vector<double> Xr(0);
  std::vector<std::array<double, 4>> X(0);

  X.reserve(10 * L);
  Xl.reserve(10 * L);
  Xr.reserve(10 * L);

  std::cout << "\nBurn-in (" << B << " generations)...";
  std::cout.flush();

  for(size_t i = 0; i < B; i++)
  {
    evolve(X, Xl, Xr, Ne, u, r, L, s, gen);
    cleanup(X, Xl, Xr);
  }

  std::cout << "done.\n";
  std::cout << "\nEvolving population(s) (" << G << " generations)...\n";

  double sum_Hl = 0.;
  double sum_Hr = 0.;
  double sum_pi2 = 0.;
  double sum_Dz = 0.;
  double sum_Dsqr = 0.;

  for(size_t i = 0; i < G; ++i)
  {
    evolve(X, Xl, Xr, Ne, L, u, r, s, gen);
    cleanup(X, Xl, Xr);

    //printX(X);
    //printXl(Xl);
    //printXr(Xr);

    sum_Hl += getSumHl(Xl);
    sum_Hr += getSumHr(Xr);
    sum_Dz += getSumDz(X);
    sum_Dsqr += getSumDsqr(X);
    sum_pi2 += getSumPi2(X);

    if(i % Ne == 1)
    {
      std::cout << "Generation " << i << "\n";
      std::cout << "\tavg_Hl = " << sum_Hl / L / i << "\n";
      std::cout << "\tavg_Hr = " << sum_Hr / L / i << "\n";
      std::cout << "\tavg_Dz = " << sum_Dz / L / i << "\n";
      std::cout << "\tavg_Dsqr = " << sum_Dsqr / L / i << "\n";
      std::cout << "\tavg_pi2 = " << sum_pi2 / L / i << "\n\n";
    }
  }

  std::cout << "avg_Hl = " << sum_Hl / L / G << "\n";
  std::cout << "avg_Hr = " << sum_Hr / L / G << "\n";
  std::cout << "avg_Dz = " << sum_Dz / L / G << "\n";
  std::cout << "avg_Dsqr = " << sum_Dsqr / L / G << "\n";
  std::cout << "avg_pi2 = " << sum_pi2 / L / G << "\n\n";

  gsl_rng_free(gen);

  twoLocusSim.done();
  return 0;
}
