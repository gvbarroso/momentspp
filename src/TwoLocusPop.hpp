/*
 * Authors: Gustavo V. Barroso
 * Created: 24/10/2023
 * Last modified: 24/10/2023
 */


#ifndef _TWO_LOCUS_POP_H_
#define _TWO_LOCUS_POP_H_

#include <algorithm>
#include <utility>
#include <vector>
#include <iostream>
#include <array>

#include "TwoLocusPair.hpp"

class TwoLocusPop
{

private:
  size_t id_;
  size_t l_; // length of segment
  unsigned int n_; // number of haplotypes
  std::vector<TwoLocusPair> pairs_;

public:
  TwoLocusPop():
  id_(0),
  l_(1),
  n_(0),
  pairs_(0)
  { }

  TwoLocusPop(size_t id, size_t l, size_t n, const std::vector<TwoLocusPair>& pairs):
  id_(id),
  l_(l),
  n_(n),
  pairs_(pairs)
  { }

public:
  void printAttributes(std::ostream& stream)
  {
    stream << "Pop. id: " << id_ << " (" << pairs_.size() << " segregating pairs)\n";

    std::array<double, 4> props = fetchAvgProps();
    std::array<double, 5> stats = fetchAvgStats();

    stream << "avg f_ab = " << props[0] << ", ";
    stream << "avg f_Ab = " << props[1] << ", ";
    stream << "avg f_aB = " << props[2] << ", ";
    stream << "avg f_AB = " << props[3] << "\n";

    stream << "avg Hl = " << stats[0] << ", ";
    stream << "avg Hr = " << stats[1] << ", ";
    stream << "avg Pi2 = " << stats[2] << ", ";
    stream << "avg Dz = " << stats[3] << ", ";
    stream << "avg D^2 = " << stats[4] << "\n\n";
  }

  unsigned int getN()
  {
    return n_;
  }

  size_t getId()
  {
    return id_;
  }

  size_t getL()
  {
    return l_;
  }

  void setId(size_t id)
  {
    id_ = id;
  }

  void setPopSize(unsigned int n)
  {
    n_ = n;

    for(auto it = std::begin(pairs_); it != std::end(pairs_); ++it)
      it->setPopSize(n);
  }

  const std::vector<TwoLocusPair>& getPairs() const
  {
    return pairs_;
  }

  std::vector<TwoLocusPair>& getPairs()
  {
    return pairs_;
  }

  std::array<double, 4> fetchAvgProps()
  {
    std::array<double, 4> props = { 0, 0, 0, 0 };

    if(pairs_.size() > 0)
    {
      for(auto it = std::begin(pairs_); it != std::end(pairs_); ++it)
      {
        props[0] += it->getProp_ab();
        props[1] += it->getProp_Ab();
        props[2] += it->getProp_aB();
        props[3] += it->getProp_AB();
      }

      props[0] /= pairs_.size();
      props[1] /= pairs_.size();
      props[2] /= pairs_.size();
      props[3] /= pairs_.size();
    }

    return props;
  }

  // within-population stats
  std::array<double, 5> fetchAvgStats()
  {
    std::array<double, 5> stats = { 0, 0, 0, 0, 0 };

    if(pairs_.size() > 0)
    {
      for(auto it = std::begin(pairs_); it != std::end(pairs_); ++it)
      {
        stats[0] += it->fetchHl();
        stats[1] += it->fetchHr();

        if(it->bothPolymorphic())
        {
          stats[2] += it->fetchPi2();
          stats[3] += it->fetchDz();
          stats[4] += it->fetchDsqr();
        }
      }

      stats[0] /= l_;
      stats[1] /= l_;
      stats[2] /= l_;
      stats[3] /= l_;
      stats[4] /= l_;
    }

    return stats;
  }

  void evolve_random(const gsl_rng* gen, size_t g, double u, double r, double s)
  {
    size_t unlinkedMuts = gsl_ran_poisson(gen, l_ * n_ * u);

    for(size_t j = 0; j < unlinkedMuts; ++j)
    {
      unsigned int c_ab = n_ - 1;
      unsigned int c_Ab = 0;
      unsigned int c_aB = 0;
      unsigned int c_AB = 0;

      if(gsl_rng_uniform(gen) < 0.5)
        c_Ab = 1;

      else
        c_aB = 1;

      TwoLocusPair newPair(g, c_ab, c_Ab, c_aB, c_AB);
      pairs_.emplace_back(newPair);
    }

    for(auto it = std::begin(pairs_); it != std::end(pairs_);)
    {
      it->evolve_random(gen, g, u, r, s);

      if(it->monomorphic() && it->mutatedBoth())
        it = pairs_.erase(it);

      else
        ++it;
    }
  }

  void evolve_random(TwoLocusPop& other, const gsl_rng* gen, size_t g, double m, double u, double r, double s)
  {
    size_t unlinkedMuts = gsl_ran_poisson(gen, l_ * n_ * u);

    for(size_t j = 0; j < unlinkedMuts; ++j)
    {
      unsigned int c_ab = n_ - 1;
      unsigned int c_Ab = 0;
      unsigned int c_aB = 0;
      unsigned int c_AB = 0;

      if(gsl_rng_uniform(gen) < 0.5)
        c_Ab = 1;

      else
        c_aB = 1;

      TwoLocusPair newPair(g, c_ab, c_Ab, c_aB, c_AB);
      pairs_.emplace_back(newPair);
    }

    migrate_random_(other, gen, m);

    for(auto it = std::begin(pairs_); it != std::end(pairs_);)
    {
      it->evolve_random(gen, g, u, r, s);

      if(it->monomorphic() && it->mutatedBoth())
        it = pairs_.erase(it);

      else
        ++it;
    }
  }

  void evolve_det(const gsl_rng* gen, size_t g, double u, double r, double s)
  {
    size_t unlinkedMuts = gsl_ran_poisson(gen, l_ * n_ * u);

    for(size_t j = 0; j < unlinkedMuts; ++j)
    {
      unsigned int c_ab = n_ - 1;
      unsigned int c_Ab = 0;
      unsigned int c_aB = 0;
      unsigned int c_AB = 0;

      if(gsl_rng_uniform(gen) < 0.5)
        c_Ab = 1;

      else
        c_aB = 1;

      TwoLocusPair newPair(g, c_ab, c_Ab, c_aB, c_AB);
      pairs_.emplace_back(newPair);
    }

    for(auto it = std::begin(pairs_); it != std::end(pairs_);)
    {
      it->evolve_det(gen, g, u, r, s);

      if(it->monomorphic() && it->mutatedBoth())
        it = pairs_.erase(it);

      else
        ++it;
    }
  }

  void evolve_det(TwoLocusPop& other, const gsl_rng* gen, size_t g, double m, double u, double r, double s)
  {
    size_t unlinkedMuts = gsl_ran_poisson(gen, l_ * n_ * u);

    for(size_t j = 0; j < unlinkedMuts; ++j)
    {
      unsigned int c_ab = n_ - 1;
      unsigned int c_Ab = 0;
      unsigned int c_aB = 0;
      unsigned int c_AB = 0;

      if(gsl_rng_uniform(gen) < 0.5)
        c_Ab = 1;

      else
        c_aB = 1;

      TwoLocusPair newPair(g, c_ab, c_Ab, c_aB, c_AB);
      pairs_.emplace_back(newPair);
    }

    migrate_det_(other, m);

    for(auto it = std::begin(pairs_); it != std::end(pairs_);)
    {
      it->evolve_det(gen, g, u, r, s);

      if(it->monomorphic() && it->mutatedBoth())
        it = pairs_.erase(it);

      else
        ++it;
    }
  }

private:
  void migrate_random_(TwoLocusPop& other, const gsl_rng* gen, double m);

  void migrate_det_(TwoLocusPop& other, double m);
};

#endif

