/*
 * Author: Gustavo V. Barroso
 * Created: 20/10/2023
 * Last modified: 24/10/2024
 * Source code for twoLocusSim
 *
 */

////////////////////////////////////////////////////////////////////////////////////////////////////

#include <vector>
#include <limits>
#include <thread>
#include <array>
#include <numeric>
#include <random>

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
  std::cout << "* Authors: G. V. Barroso                 Last Modif. 23/Oct/2024 *" << std::endl;
  std::cout << "*          A. P. Ragsdale                                        *" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "******************************************************************" << std::endl;

  std::cout << "\nCompiled on: " << __DATE__ << std::endl;
  std::cout << "Compiled at: " << __TIME__ << std::endl << std::endl;

  if(argc == 1)
  {
    std::cout << "To use TwoLocusSim, fill in a text file with the following options and execute from the command line:\ntwolocussim params=file_name\n\n";

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
  double u = bpp::ApplicationTools::getParameter<double>("u", params, 1e-8, "", 0);
  double r = bpp::ApplicationTools::getParameter<double>("r", params, 1e-7, "", 0);
  double s = bpp::ApplicationTools::getParameter<double>("s", params, 0., "", 0);
  std::string label = bpp::ApplicationTools::getStringParameter("label", params, "test", "", true, 4);

  std::cout << "\nSimulation setup:\n";
  std::cout << "\t" << L << " loci\n";
  std::cout << "\tfor " << G << " generations\n";
  std::cout << "\tNe = " << Ne << "\n";
  std::cout << "\tu = " << u << "\n";
  std::cout << "\tr = " << r << "\n";
  std::cout << "\ts = " << s << "\n";

  unsigned int B = 40 * Ne;

  // helper generator to avoid setting same seed when replicates start at the "same time"
  std::mt19937 rng;
  std::array<int, std::mt19937::state_size> seedData;
  unsigned sem = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine re(sem);
  std::generate_n(seedData.data(), seedData.size(), std::ref(re));
  std::seed_seq seq(std::begin(seedData), std::end(seedData));
  rng.seed(seq);

  std::poisson_distribution<> pois(100000000);

  struct timeval tv;
  gettimeofday(&tv, 0);
  unsigned long seed = tv.tv_sec + tv.tv_usec + pois(rng);

  std::cout << "random seed = " << seed << "\n";

  const gsl_rng_type * T;
  gsl_rng* gen;

  gsl_rng_env_setup();

  T = gsl_rng_default;
  gen = gsl_rng_alloc(T);

  gsl_rng_set(gen, seed);

  /*TwoLocusPop root(0, L, Ne);

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
  std::cout << "Evolving population(s) (" << G << " generations)..."; std::cout.flush();

  boost::iostreams::filtering_ostream stream_D;

  std::string tbl_file = "tbl_D_" + label + "_" + bpp::TextTools::toString(seed) + ".txt";
  std::ofstream file_D;
  file_D.open(tbl_file + ".gz", std::ios_base::out | std::ios_base::binary);

  stream_D.push(boost::iostreams::gzip_compressor());
  stream_D.push(file_D);
  stream_D << "D\tp\tq\n";


  double sum_Hl = 0.;
  double sum_Hr = 0.;
  double sum_pi2 = 0.;
  double sum_D = 0.;
  double sum_Dr = 0.;
  double sum_Dl = 0.;
  double sum_Dz = 0.;
  double sum_Dsqr = 0.;*/

  boost::iostreams::filtering_ostream stream_hap;

  std::string hap_file = "hap_trajectories_" + label + "_" + bpp::TextTools::toString(seed) + ".txt";
  std::ofstream file_hap;
  file_hap.open(hap_file + ".gz", std::ios_base::out | std::ios_base::binary);

  stream_hap.push(boost::iostreams::gzip_compressor());
  stream_hap.push(file_hap);

  TwoLocusPop p1(0, L, Ne);

  for(size_t i = 0; i < G; ++i)
  {
    p1.evolve(u, r, s, gen);
    p1.cleanup(stream_hap);

    if(p1.getX().size() > 0)
      p1.printX(stream_hap);

    /*p1.computeStats();

    sum_Hl += p1.getSumHl();
    sum_Hr += p1.getSumHr();
    sum_D += p1.getSumD();
    sum_Dr += p1.getSumDr();
    sum_Dl += p1.getSumDl();
    sum_Dz += p1.getSumDz();
    sum_Dsqr += p1.getSumDsqr();
    sum_pi2 += p1.getSumPi2();

    if(i % Ne == 0)
      p1.tabulate_Ds(stream_D);*/

  }

  //boost::iostreams::close(stream_D);
  boost::iostreams::close(stream_hap);

  /*std::string file = "stats_" + label + "_" + bpp::TextTools::toString(seed) + ".txt";
  std::ofstream fout(file);

  fout << "Hl = " << sum_Hl / L / G << "\n";
  fout << "Hr = " << sum_Hr / L / G << "\n";
  fout << "D = " << sum_D / L / L / G << "\n";
  fout << "Dr = " << sum_Dr / L / L / G << "\n";
  fout << "Dl = " << sum_Dl / L / L / G << "\n";
  fout << "Dz = " << sum_Dz / L / L / G << "\n";
  fout << "D^2 = " << sum_Dsqr / L / L / G << "\n";
  fout << "pi2 = " << sum_pi2 / L / L / G << "\n";
  fout << "random_seed = " << seed << "\n";

  fout.close();*/

  //std::cout << "done.\nCheck output file " << file << ".\n\n";

  gsl_rng_free(gen);

  twoLocusSim.done();
  return 0;
}
