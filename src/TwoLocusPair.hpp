/*
 * Authors: Gustavo V. Barroso
 * Created: 19/10/2023
 * Last modified: 08/11/2023
 */


#ifndef _TWO_LOCUS_PAIR_H_
#define _TWO_LOCUS_PAIR_H_

#include <algorithm>
#include <utility>
#include <cassert>
#include <chrono>
#include <iostream>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include <Bpp/Exceptions.h>

class TwoLocusPair
{

private:
  // counts don't have to be integers here (for deterministic sel and rec)
  double count_ab_;
  double count_Ab_;
  double count_aB_;
  double count_AB_;

  bool mutatedLeft_;
  bool mutatedRight_;

public:
  TwoLocusPair():
  count_ab_(0),
  count_Ab_(0),
  count_aB_(0),
  count_AB_(0),
  mutatedLeft_(0),
  mutatedRight_(0)
  { }

  TwoLocusPair(double count_ab, double count_Ab, double count_aB, double count_AB):
  count_ab_(count_ab),
  count_Ab_(count_Ab),
  count_aB_(count_aB),
  count_AB_(count_AB),
  mutatedLeft_(0.),
  mutatedRight_(0.)
  {
    assert(!mutatedBoth());

    mutatedLeft_ = fetchP() > 0.;
    mutatedRight_ = fetchQ() > 0.;
  }

public:
  void printAttributes(std::ostream& stream)
  {
    stream << std::fixed << std::showpoint;
    stream << std::setprecision(6);
    stream << "c_ab = " << count_ab_ << ", ";
    stream << "c_Ab = " << count_Ab_ << ", ";
    stream << "c_aB = " << count_aB_ << ", ";
    stream << "c_AB = " << count_AB_ << "\n";
  }

  double getCount_ab()
  {
    return count_ab_;
  }

  double getCount_Ab()
  {
    return count_Ab_;
  }

  double getCount_aB()
  {
    return count_aB_;
  }

  double getCount_AB()
  {
    return count_AB_;
  }

  bool mutatedLeft()
  {
    return mutatedLeft_;
  }

  bool mutatedRight()
  {
    return mutatedRight_;
  }

  bool mutatedBoth()
  {
    return mutatedLeft_ && mutatedRight_;
  }

  bool bothPolymorphic()
  {
    return (fetchP() > 0. && fetchP() < 1. && fetchQ() > 0. && fetchQ() < 1.);
  }

  bool monomorphic()
  {
    double n = getNumHaps();
    return (count_ab_ == n || count_Ab_ == n || count_aB_ == n || count_AB_ == n);
  }

  double getNumHaps()
  {
    return count_ab_ + count_Ab_ + count_aB_ + count_AB_;
  }

  double fetchP()
  {
    double n = getNumHaps();
    return (count_Ab_ + count_AB_) / n;
  }

  double fetchQ()
  {
    double n = getNumHaps();
    return (count_aB_ + count_AB_) / n;
  }

  double fetchD()
  {
    double n = getNumHaps();
    return (count_ab_ * count_AB_ - count_Ab_ * count_aB_) / (n * n);
  }

  double fetchHl()
  {
    return fetchP() * (1. - fetchP());
  }

  double fetchHr()
  {
    return fetchQ() * (1. - fetchQ());
  }

  double fetchDz()
  {
    return fetchD() * (1. - 2. * fetchP()) * (1. - 2. * fetchQ());
  }

  double fetchDsqr()
  {
    double n = getNumHaps();
    return ((count_ab_ * count_AB_ - count_Ab_ * count_aB_) * (count_ab_ * count_AB_ - count_Ab_ * count_aB_)) / (n * n * n * n);
  }

  double fetchPi2()
  {
    return fetchHl() * fetchHr();
  }

  void setPopSize(double n)
  {
    count_ab_ *= n;
    count_Ab_ *= n;
    count_aB_ *= n;
    count_AB_ *= n;
  }

  void evolve_random(const gsl_rng* gen, double u, double r, double s)
  {
    if(!mutatedBoth())
      mutate_(gen, u);

    recombineRandom_(gen, r);
    selectRandom_(gen, s);
    drift_(gen);
  }

  void evolve_det(const gsl_rng* gen, double u, double r, double s)
  {
    if(!mutatedBoth())
      mutate_(gen, u);

    recombineDet_(r);
    selectDet_(s);
    drift_(gen);
  }

private:
  // to mutate monomorphic (unmutated) locus (either left or right)
  void mutate_(const gsl_rng* gen, double u)
  {
    assert(mutatedBoth() == false);

    double n_haps = getNumHaps();

    if(gsl_rng_uniform(gen) < n_haps * u)
      mutateOtherLocus_(gen);
  }

void mutateOtherLocus_(const gsl_rng* gen)
  {
    assert(mutatedBoth() == false);

    double n_haps = getNumHaps();

    if(mutatedLeft_)
    {
      if(gsl_rng_uniform(gen) < count_Ab_ / n_haps)
      {
        --count_Ab_;
        ++count_AB_;
      }

      else
      {
        --count_ab_;
        ++count_aB_;
      }

      mutatedRight_ = true;
    }

    else if(mutatedRight_)
    {
      if(gsl_rng_uniform(gen) < count_aB_ / n_haps)
      {
        --count_aB_;
        ++count_AB_;
      }

      else
      {
        --count_ab_;
        ++count_Ab_;
      }

      mutatedLeft_ = true;
    }
  }

  void drift_(const gsl_rng* gen)
  {
    unsigned int n = std::round(getNumHaps());
    //std::cout << "N before drift: " << std::setprecision(12) << n << "\n";
    //printAttributes(std::cout);
    unsigned int next[4];
    double probs[4] = { count_ab_, count_Ab_, count_aB_, count_AB_ };

    gsl_ran_multinomial(gen, 4, n, probs, next);

    count_ab_ = next[0];
    count_Ab_ = next[1];
    count_aB_ = next[2];
    count_AB_ = next[3];
    //std::cout << "N after drift: " << std::setprecision(12) << getNumHaps() << "\n";
    //printAttributes(std::cout);
  }

  void recombineRandom_(const gsl_rng* gen, double r)
  {
    double p = fetchP();
    double q = fetchQ();

    unsigned int count_ab_rec = gsl_ran_binomial(gen, r, count_ab_);
    unsigned int count_Ab_rec = gsl_ran_binomial(gen, r, count_Ab_);
    unsigned int count_aB_rec = gsl_ran_binomial(gen, r, count_aB_);
    unsigned int count_AB_rec = gsl_ran_binomial(gen, r, count_AB_);

    unsigned int ab_to_Ab = gsl_ran_binomial(gen, p, count_ab_rec);
    ab_to_Ab = gsl_ran_binomial(gen, 0.5, ab_to_Ab);

    unsigned int ab_to_aB = gsl_ran_binomial(gen, q, count_ab_rec);
    ab_to_aB = gsl_ran_binomial(gen, 0.5, ab_to_aB);

    unsigned int ab_stay_ab = count_ab_ - ab_to_Ab - ab_to_aB;

    unsigned int Ab_to_ab = gsl_ran_binomial(gen, 1. - p, count_Ab_rec);
    Ab_to_ab = gsl_ran_binomial(gen, 0.5, Ab_to_ab);

    unsigned int Ab_to_AB = gsl_ran_binomial(gen, q, count_Ab_rec);
    Ab_to_AB = gsl_ran_binomial(gen, 0.5, Ab_to_AB);

    unsigned int Ab_stay_Ab = count_Ab_ - Ab_to_ab - Ab_to_AB;

    unsigned int aB_to_AB = gsl_ran_binomial(gen, p, count_aB_rec);
    aB_to_AB = gsl_ran_binomial(gen, 0.5, aB_to_AB);

    unsigned int aB_to_ab = gsl_ran_binomial(gen, 1. - q, count_aB_rec);
    aB_to_ab = gsl_ran_binomial(gen, 0.5, aB_to_ab);

    unsigned int aB_stay_aB = count_aB_ - aB_to_AB - aB_to_ab;

    unsigned int AB_to_aB = gsl_ran_binomial(gen, 1. - p, count_AB_rec);
    AB_to_aB = gsl_ran_binomial(gen, 0.5, AB_to_aB);

    unsigned int AB_to_Ab = gsl_ran_binomial(gen, 1. - q, count_AB_rec);
    AB_to_Ab = gsl_ran_binomial(gen, 0.5, AB_to_Ab);

    unsigned int AB_stay_AB = count_AB_ - AB_to_aB - AB_to_Ab;

    count_ab_ = ab_stay_ab + Ab_to_ab + aB_to_ab;
    count_Ab_ = Ab_stay_Ab + ab_to_Ab + AB_to_Ab;
    count_aB_ = aB_stay_aB + ab_to_aB + AB_to_aB;
    count_AB_ = AB_stay_AB + Ab_to_AB + aB_to_AB;
  }

  void selectRandom_(const gsl_rng* gen, double s)
  {
    double p[4] = { count_ab_, count_Ab_, count_aB_, count_AB_ };
    unsigned int Ab_replacement[4] = { 0, 0, 0, 0 };
    unsigned int AB_replacement[4] = { 0, 0, 0, 0 };

    unsigned int count_Ab_dead = gsl_ran_binomial(gen, std::abs(s), count_Ab_);
    unsigned int count_AB_dead = gsl_ran_binomial(gen, std::abs(s), count_AB_);

    if(count_Ab_dead > 0)
      gsl_ran_multinomial(gen, 4, count_Ab_dead, p, Ab_replacement);

    if(count_AB_dead > 0)
      gsl_ran_multinomial(gen, 4, count_AB_dead, p, AB_replacement);

    count_ab_ = count_ab_ + Ab_replacement[0] + AB_replacement[0];
    count_Ab_ = count_Ab_ - count_Ab_dead + Ab_replacement[1] + AB_replacement[1];
    count_aB_ = count_aB_ + Ab_replacement[2] + AB_replacement[2];
    count_AB_ = count_AB_ - count_AB_dead + Ab_replacement[3] + AB_replacement[3];
  }

  void recombineDet_(double r)
  {
    //std::cout << "N before rec: " << std::setprecision(12) << getNumHaps() << "\n";
    //printAttributes(std::cout);
    count_ab_ -= r * fetchD();
    count_Ab_ += r * fetchD();
    count_aB_ += r * fetchD();
    count_AB_ -= r * fetchD();
    //std::cout << "N after rec: " << std::setprecision(12) << getNumHaps() << "\n";
    //printAttributes(std::cout);
  }

  void selectDet_(double s)
  {
    double n = getNumHaps();
    //std::cout << "N before sel: " << std::setprecision(12) << n << "\n";
    //printAttributes(std::cout);
    count_ab_ *= (1. - s * fetchP());
    count_Ab_ *= (1. + s * (1. - fetchP()));
    count_aB_ *= (1. - s * fetchP());
    count_AB_ *= (1. + s * (1. - fetchP()));

    double f = getNumHaps();

    //std::cout << "N after  sel: " << std::setprecision(12) << f << "\n";

    count_ab_ *= n / f;
    count_Ab_ *= n / f;
    count_aB_ *= n / f;
    count_AB_ *= n / f;

    //std::cout << "N after norm: " << std::setprecision(12) << getNumHaps() << "\n";
    //printAttributes(std::cout);
  }
};

#endif

