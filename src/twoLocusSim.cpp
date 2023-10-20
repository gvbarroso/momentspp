/*
 * Author: Gustavo V. Barroso
 * Created: 20/10/2023
 * Last modified: 20/10/2023
 * Source code for twoLocusSim
 *
 */

////////////////////////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <vector>
#include <map>
#include <limits>
#include <thread>

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
  std::cout << "* Authors: G. V. Barroso                 Last Modif. 20/Oct/2023 *" << std::endl;
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

  size_t L = bpp::ApplicationTools::getParameter<size_t>("L", params, 1e+4, "", 0);
  size_t N = bpp::ApplicationTools::getParameter<size_t>("N", params, 1e+4, "", 0);
  double u = bpp::ApplicationTools::getParameter<double>("u", params, 1e-6, "", 0);
  double r = bpp::ApplicationTools::getParameter<double>("r", params, 1e-7, "", 0);
  double s = bpp::ApplicationTools::getParameter<double>("s", params, -1e-4, "", 0);

  std::array<int, 624> seedData;
  unsigned sem = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine re(sem);
  std::generate_n(seedData.data(), seedData.size(), std::ref(re));
  std::seed_seq seq(std::begin(seedData), std::end(seedData));

  std::mt19937 gen;
  gen.seed(seq);

  std::uniform_real_distribution<double> unif(0., 1.);
  std::poisson_distribution<size_t> pois(L * N * u);
  std::vector<TwoLocusPair> pairs(0);
  pairs.reserve(L);

  // TODO include vector of N's and vector of G's as options to set discrete epochs
  size_t numGen = 1e+6;
  size_t mutablePairs = 0;

  for(size_t i = 0; i < numGen; ++i) // for each epoch, from past to present
  {
    size_t mutCount = pois(gen);
    size_t mutsInPairs = mutCount * (mutablePairs / L);
    size_t unlinkedMuts = mutCount - mutsInPairs;

    std::shuffle(std::begin(pairs), std::end(pairs), gen);

    size_t track = 0;
    for(size_t j = 0; j < pairs.size(); ++j)
    {
      if(track < mutsInPairs)
      {
        if(!pairs[j].mutatedBoth())
        {
          pairs[j].mutate();
          ++track;
        }
      }

      else
        break;
    }

    for(size_t j = 0; j < unlinkedMuts; ++j)
    {
      size_t c_ab = N - 1;
      size_t c_Ab = 0;
      size_t c_aB = 0;
      size_t c_AB = 0;

      if(unif(gen) < 0.5)
        c_Ab = 1;

      else
        c_aB = 1;

      TwoLocusPair newPair(c_ab, c_Ab, c_aB, c_AB);
      pairs.emplace_back(newPair);
    }
  }



  twoLocusSim.done();
  return 0;
}
