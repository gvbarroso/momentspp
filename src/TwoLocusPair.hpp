/*
 * Authors: Gustavo V. Barroso
 * Created: 19/10/2023
 * Last modified: 24/10/2023
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
  size_t gen_origin_;

  // for random operators
  unsigned int count_ab_;
  unsigned int count_Ab_;
  unsigned int count_aB_;
  unsigned int count_AB_;

  unsigned int n_; // number of haplotypes

  // for deterministic operators
  double prop_ab_;
  double prop_Ab_;
  double prop_aB_;
  double prop_AB_;

  bool mutatedLeft_;
  bool mutatedRight_;

public:
  TwoLocusPair():
  gen_origin_(0),
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

  TwoLocusPair(size_t gen_origin, unsigned int count_ab, unsigned int count_Ab, unsigned int count_aB, unsigned int count_AB):
  gen_origin_(gen_origin),
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
  }

public:
  void printAttributes(std::ostream& stream)
  {
    stream << "c_ab = " << count_ab_ << " (" << prop_ab_ << "), ";
    stream << "c_Ab = " << count_Ab_ << " (" << prop_Ab_ << "), ";
    stream << "c_aB = " << count_aB_ << " (" << prop_aB_ << "), ";
    stream << "c_AB = " << count_AB_ << " (" << prop_AB_ << ")\n";
  }

  size_t getGenOrigin()
  {
    return gen_origin_;
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

  bool bothPolymorphic()
  {
    return (fetchP() > 0. && fetchP() < 1.) && (fetchQ() > 0. && fetchQ() < 1.);
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

  void evolve_random(const gsl_rng* gen, double u, double r, double s)
  {
    drift_(gen);

    if(!mutatedBoth())
      mutate_(gen, u);

    recombineRandom_(gen, r);
    selectRandom_(gen, s);
  }

  void evolve_det(const gsl_rng* gen, double u, double r, double s)
  {
    drift_(gen);

    if(!mutatedBoth())
      mutate_(gen, u);

    recombineDet_(r);
    selectDet_(s);
  }

private:
  void updateProps_()
  {
    prop_ab_ = static_cast<double>(count_ab_) / static_cast<double>(n_);
    prop_Ab_ = static_cast<double>(count_Ab_) / static_cast<double>(n_);
    prop_aB_ = static_cast<double>(count_aB_) / static_cast<double>(n_);
    prop_AB_ = static_cast<double>(count_AB_) / static_cast<double>(n_);

    normalize_();
  }

  void normalize_()
  {
    double sum = prop_ab_ + prop_Ab_ + prop_aB_ + prop_AB_;

    prop_ab_ /= sum;
    prop_Ab_ /= sum;
    prop_aB_ /= sum;
    prop_AB_ /= sum;
  }

  // used only once, to mutate monomorphic locus (either left or right)
  void mutate_(const gsl_rng* gen, double u)
  {
    assert(mutatedBoth() == false);

    if(gsl_rng_uniform(gen) < n_ * u)
    {
      if(mutatedLeft_)
      {
        if(gsl_rng_uniform(gen) < prop_Ab_)
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
        if(gsl_rng_uniform(gen) < prop_aB_)
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

    updateProps_();
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

    updateProps_();
  }

  void recombineRandom_(const gsl_rng* gen, double r)
  {
    assert(count_ab_ + count_Ab_ + count_aB_ + count_AB_ == n_);

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

    unsigned int aB_to_AB = gsl_ran_binomial(gen, q, count_aB_rec);
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

    updateProps_();
  }

  void selectRandom_(const gsl_rng* gen, double s)
  {
    assert(count_ab_ + count_Ab_ + count_aB_ + count_AB_ == n_);

    double p[4] = { prop_ab_, prop_Ab_, prop_aB_, prop_AB_ };
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

    assert(count_ab_ + count_Ab_ + count_aB_ + count_AB_ == n_);

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

