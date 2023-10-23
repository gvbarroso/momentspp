/*
 * Authors: Gustavo V. Barroso
 * Created: 19/10/2023
 * Last modified: 23/10/2023
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

class TwoLocusPair
{

private:
  // for random operators
  unsigned int count_ab_;
  unsigned int count_Ab_;
  unsigned int count_aB_;
  unsigned int count_AB_;

  unsigned int n_;

  // for deterministic operators
  double prop_ab_;
  double prop_Ab_;
  double prop_aB_;
  double prop_AB_;

  bool mutatedLeft_;
  bool mutatedRight_;

public:
  TwoLocusPair():
  count_ab_(0),
  count_Ab_(0),
  count_aB_(0),
  count_AB_(0),
  n_(0),
  prop_ab_(0),
  prop_Ab_(0),
  prop_aB_(0),
  prop_AB_(0),
  mutatedLeft_(0),
  mutatedRight_(0)
  { }

  TwoLocusPair(unsigned int count_ab, unsigned int count_Ab, unsigned int count_aB, unsigned int count_AB):
  count_ab_(count_ab),
  count_Ab_(count_Ab),
  count_aB_(count_aB),
  count_AB_(count_AB),
  n_(count_ab + count_Ab + count_aB + count_AB),
  prop_ab_(static_cast<double>(count_ab) / n_),
  prop_Ab_(static_cast<double>(count_Ab) / n_),
  prop_aB_(static_cast<double>(count_aB) / n_),
  prop_AB_(static_cast<double>(count_AB) / n_),
  mutatedLeft_(fetchP() > 0.),
  mutatedRight_(fetchQ() > 0.)
  {
    assert(!mutatedBoth());

    printAttributes(std::cout);
  }

public:
  void printAttributes(std::ostream& stream)
  {
    stream << "c_ab = " << count_ab_ << " (" << prop_ab_ << ")\n";
    stream << "c_Ab = " << count_Ab_ << " (" << prop_Ab_ << ")\n";
    stream << "c_aB = " << count_aB_ << " (" << prop_aB_ << ")\n";
    stream << "c_AB = " << count_AB_ << " (" << prop_AB_ << ")\n";
  }

  unsigned int getCount_ab()
  {
    return count_ab_;
  }

  unsigned int getCount_Ab()
  {
    return count_Ab_;
  }

  unsigned int getCount_aB()
  {
    return count_aB_;
  }

  unsigned int getCount_AB()
  {
    return count_AB_;
  }

  unsigned int getProp_ab()
  {
    return prop_ab_;
  }

  unsigned int getProp_Ab()
  {
    return prop_Ab_;
  }

  unsigned int getProp_aB()
  {
    return prop_aB_;
  }

  unsigned int getProp_AB()
  {
    return prop_AB_;
  }

  unsigned int getN()
  {
    return n_;
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

  bool monomorphic()
  {
    return (count_ab_ == n_ || count_Ab_ == n_ || count_aB_ == n_ || count_AB_ == n_);
  }

  double fetchP()
  {
    return prop_Ab_ + prop_AB_;
  }

  double fetchQ()
  {
    return prop_aB_ + prop_AB_;
  }

  double fetchD()
  {
    return prop_ab_ * prop_AB_ - prop_Ab_ * prop_aB_;
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
    return (prop_ab_ * prop_AB_ - prop_Ab_ * prop_aB_) * (prop_ab_ * prop_AB_ - prop_Ab_ * prop_aB_);
  }

  double fetchPi2()
  {
    return fetchP() * (1. - fetchP()) * fetchQ() * (1. - fetchQ());
  }

  void evolve_random(const gsl_rng* gen, double r, double s)
  {
    drift_(gen);
    recombineRandom_(gen, r);
    selectRandom_(gen, s);
  }

  // used only once, to mutate monomorphic locus (either left or right)
  void mutate(const gsl_rng* gen)
  {
    assert(mutatedBoth() == false);

    if(mutatedLeft_)
    {
      if(gsl_rng_uniform(gen) < static_cast<double>(count_Ab_) / n_)
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
      if(gsl_rng_uniform(gen) < static_cast<double>(count_aB_) / n_)
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

private:
  void updateProps_()
  {
    prop_ab_ = static_cast<double>(count_ab_) / static_cast<double>(n_);
    prop_Ab_ = static_cast<double>(count_Ab_) / static_cast<double>(n_);
    prop_aB_ = static_cast<double>(count_aB_) / static_cast<double>(n_);
    prop_AB_ = static_cast<double>(count_AB_) / static_cast<double>(n_);
  }

  void normalize_()
  {
    double sum = prop_ab_ + prop_Ab_ + prop_aB_ + prop_AB_;

    prop_ab_ /= sum;
    prop_Ab_ /= sum;
    prop_aB_ /= sum;
    prop_AB_ /= sum;
  }

  void drift_(const gsl_rng* gen)
  {
    assert(count_ab_ + count_Ab_ + count_aB_ + count_AB_ == n_);

    unsigned int next[4] = { count_ab_, count_Ab_, count_aB_, count_AB_ };
    double probs[4] = { prop_ab_, prop_Ab_, prop_aB_, prop_AB_ };

    gsl_ran_multinomial(gen, 4, n_, probs, next);

    count_ab_ = next[0];
    count_Ab_ = next[1];
    count_aB_ = next[2];
    count_AB_ = next[3];

    printAttributes(std::cout);
    updateProps_();
  }

  void recombineRandom_(const gsl_rng* gen, double r)
  {
    assert(count_ab_ + count_Ab_ + count_aB_ + count_AB_ == n_);

    double p = fetchP();
    double q = fetchQ();

    unsigned int count_ab_rec = gsl_ran_binomial(gen, count_ab_, r);
    unsigned int count_Ab_rec = gsl_ran_binomial(gen, count_Ab_, r);
    unsigned int count_aB_rec = gsl_ran_binomial(gen, count_aB_, r);
    unsigned int count_AB_rec = gsl_ran_binomial(gen, count_AB_, r);

    unsigned int ab_to_Ab = gsl_ran_binomial(gen, count_ab_rec, p);
    ab_to_Ab = gsl_ran_binomial(gen, ab_to_Ab, 0.5);

    unsigned int ab_to_aB = gsl_ran_binomial(gen, count_ab_rec, q);
    ab_to_aB = gsl_ran_binomial(gen, ab_to_aB, 0.5);

    unsigned int ab_stay_ab = count_ab_ - ab_to_Ab - ab_to_aB;

    unsigned int Ab_to_ab = gsl_ran_binomial(gen, count_Ab_rec, 1. - p);
    Ab_to_ab = gsl_ran_binomial(gen, Ab_to_ab, 0.5);

    unsigned int Ab_to_AB = gsl_ran_binomial(gen, count_Ab_rec, q);
    Ab_to_AB = gsl_ran_binomial(gen, Ab_to_AB, 0.5);

    unsigned int Ab_stay_Ab = count_Ab_ - Ab_to_ab - Ab_to_AB;

    unsigned int aB_to_AB = gsl_ran_binomial(gen, count_aB_rec, q);
    aB_to_AB = gsl_ran_binomial(gen, aB_to_AB, 0.5);

    unsigned int aB_to_ab = gsl_ran_binomial(gen, count_aB_rec, 1. - q);
    aB_to_ab = gsl_ran_binomial(gen, aB_to_ab, 0.5);

    unsigned int aB_stay_aB = count_aB_ - aB_to_AB - aB_to_ab;

    unsigned int AB_to_aB = gsl_ran_binomial(gen, count_AB_rec, 1. - p);
    AB_to_aB = gsl_ran_binomial(gen, AB_to_aB, 0.5);

    unsigned int AB_to_Ab = gsl_ran_binomial(gen, count_AB_rec, 1. - 1);
    AB_to_Ab = gsl_ran_binomial(gen, AB_to_Ab, 0.5);

    unsigned int AB_stay_AB = count_AB_ - AB_to_aB - AB_to_Ab;

    count_ab_ = ab_stay_ab + Ab_to_ab + aB_to_ab;
    count_Ab_ = Ab_stay_Ab + ab_to_Ab + AB_to_Ab;
    count_aB_ = aB_stay_aB + ab_to_aB + AB_to_aB;
    count_AB_ = AB_stay_AB + Ab_to_AB + aB_to_AB;

    updateProps_();
  }

  void selectRandom_(const gsl_rng* gen, double s)
  {
    assert(count_ab_ + count_Ab_ + count_aB_ + count_AB_ == n_);

    double p[4] = { prop_ab_, prop_Ab_, prop_aB_, prop_AB_ };
    unsigned int Ab_replacement[4] = { 0, 0, 0, 0 };
    unsigned int AB_replacement[4] = { 0, 0, 0, 0 };

    unsigned int count_Ab_dead = gsl_ran_binomial(gen, count_Ab_, std::abs(s));
    unsigned int count_AB_dead = gsl_ran_binomial(gen, count_AB_, std::abs(s));

    if(count_Ab_dead > 0)
      gsl_ran_multinomial(gen, n_, count_Ab_dead, p, Ab_replacement);

    if(count_AB_dead > 0)
      gsl_ran_multinomial(gen, n_, count_AB_dead, p, AB_replacement);

    count_ab_ = count_ab_ + Ab_replacement[0] + AB_replacement[0];
    count_Ab_ = count_Ab_ - count_Ab_dead + Ab_replacement[1] + AB_replacement[1];
    count_aB_ = count_aB_ + Ab_replacement[2] + AB_replacement[2];
    count_AB_ = count_AB_ - count_AB_dead + AB_replacement[3] + AB_replacement[3];

    updateProps_();
  }

  void recombineDet_(double r)
  {
    prop_ab_ = prop_ab_ - r * fetchD();
    prop_Ab_ = prop_Ab_ + r * fetchD();
    prop_aB_ = prop_aB_ + r * fetchD();
    prop_AB_ = prop_AB_ - r * fetchD();

    normalize_();
  }

  void selectDet_(double s)
  {
    prop_ab_ = prop_ab_ - s * prop_ab_ * fetchP();
    prop_Ab_ = prop_Ab_ + s * (1. - fetchP());
    prop_aB_ = prop_aB_ - s * prop_aB_ * fetchP();
    prop_AB_ = prop_AB_ + s * (1. - fetchP());

    normalize_();
  }
};

#endif

