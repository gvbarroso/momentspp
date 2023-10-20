/*
 * Authors: Gustavo V. Barroso
 * Created: 19/10/2023
 * Last modified: 20/10/2023
 */


#ifndef _TWO_LOCUS_PAIR_H_
#define _TWO_LOCUS_PAIR_H_

#include <random>
#include <algorithm>
#include <utility>
#include <cassert>
#include <array>
#include <chrono>

#include <gsl/gsl_randist.h>

class TwoLocusPair
{

private:
  // for random operators
  size_t count_ab_;
  size_t count_Ab_;
  size_t count_aB_;
  size_t count_AB_;

  size_t n_;

  // for deterministic operators
  double prop_ab_;
  double prop_Ab_;
  double prop_aB_;
  double prop_AB_;

  bool mutatedLeft_;
  bool mutatedRight_;

  std::binomial_distribution<size_t> gaussianJitter_ { 0.0, 1e-4};

public:
  TwoLocusPair(size_t count_ab, size_t count_Ab, size_t count_aB, size_t count_AB):
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
  { }

public:
  size_t getCount_ab()
  {
    return count_ab_;
  }

  size_t getCount_Ab()
  {
    return count_Ab_;
  }

  size_t getCount_aB()
  {
    return count_aB_;
  }

  size_t getCount_AB()
  {
    return count_AB_;
  }

  size_t getProp_ab()
  {
    return prop_ab_;
  }

  size_t getProp_Ab()
  {
    return prop_Ab_;
  }

  size_t getProp_aB()
  {
    return prop_aB_;
  }

  size_t getProp_AB()
  {
    return prop_AB_;
  }

  size_t getN()
  {
    return n_;
  }

  bool mutatedLeft()
  {
    return mutatedLeft_;
  }

  bool mutatedRight()
  {
    mutatedRight_;
  }

  bool mutatedBoth()
  {
    mutatedLeft_ && mutatedRight_;
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

  void evolve_random(std::mt19937& gen, double r, double s, double u)
  {
    if(!mutatedBoth())
    {

    }

    else
    {
    }

    updateProps_();
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

  void driftProps_(std::mt19937& gen)
  {
    unsigned int next[4] = { prop_ab_, prop_Ab_, prop_aB_, prop_AB_ };
    gsl_ran_multinomial(gen, n_, n_, next, next);

    prop_ab_ = next[0];
    prop_Ab_ = next[1];
    prop_aB_ = next[2];
    prop_AB_ = next[3];

    delete [] next;

    normalize_();
  }

  void driftCounts_(std::mt19937& gen)
  {
    assert(count_ab_ + count_Ab_ + count_aB_ + count_AB_ == n_);

    unsigned int next[4] = { count_ab_, count_Ab_, count_aB_, count_AB_ };
    gsl_ran_multinomial(gen, n_, n_, next, next);

    count_ab_ = next[0];
    count_Ab_ = next[1];
    count_aB_ = next[2];
    count_AB_ = next[3];

    delete [] next;
  }

  // used only once, to mutate monomorphic locus (either left or right)
  void mutate_(std::mt19937& gen, const std::uniform_real_distribution<double>& unif)
  {
    assert(mutatedBoth() == false);

    if(mutatedLeft_)
    {
      if(unif(gen) < static_cast<double>(count_Ab_) / n_)
      {
        --count_Ab_;
        ++count_AB_;
      }

      else
      {
        --count_ab_;
        ++count_aB_;
      }

      mutatedRight_ == true;
    }

    else if(mutatedRight_)
    {
      if(unif(gen) < static_cast<double>(count_aB_) / n_)
      {
        --count_aB_;
        ++count_AB_;
      }

      else
      {
        --count_ab_;
        ++count_Ab_;
      }

      mutatedLeft_ == true;
    }
  }

  void recombineRandom_(std::mt19937& gen, double r)
  {
    assert(count_ab_ + count_Ab_ + count_aB_ + count_AB_ == n_);

    double p = fetchP();
    double q = fetchQ();

    std::binomial_distribution<size_t> dab(count_ab_, r);
    std::binomial_distribution<size_t> dAb(count_Ab_, r);
    std::binomial_distribution<size_t> daB(count_aB_, r);
    std::binomial_distribution<size_t> dAB(count_AB_, r);

    size_t count_ab_rec = dab(gen);
    size_t count_Ab_rec = dAb(gen);
    size_t count_aB_rec = daB(gen);
    size_t count_AB_rec = dAB(gen);

    std::binomial_distribution<size_t> dab_rec_p(count_ab_rec, p);
    std::binomial_distribution<size_t> dAb_rec_p(count_Ab_rec, p);
    std::binomial_distribution<size_t> daB_rec_p(count_aB_rec, p);
    std::binomial_distribution<size_t> dAB_rec_p(count_AB_rec, p);

    std::binomial_distribution<size_t> dab_rec_ap(count_ab_rec, 1. - p);
    std::binomial_distribution<size_t> dAb_rec_ap(count_Ab_rec, 1. - p);
    std::binomial_distribution<size_t> daB_rec_ap(count_aB_rec, 1. - p);
    std::binomial_distribution<size_t> dAB_rec_ap(count_AB_rec, 1. - p);

    std::binomial_distribution<size_t> dab_rec_q(count_ab_rec, q);
    std::binomial_distribution<size_t> dAb_rec_q(count_Ab_rec, q);
    std::binomial_distribution<size_t> daB_rec_q(count_aB_rec, q);
    std::binomial_distribution<size_t> dAB_rec_q(count_AB_rec, q);

    std::binomial_distribution<size_t> dab_rec_aq(count_ab_rec, 1. - q);
    std::binomial_distribution<size_t> dAb_rec_aq(count_Ab_rec, 1. - q);
    std::binomial_distribution<size_t> daB_rec_aq(count_aB_rec, 1. - q);
    std::binomial_distribution<size_t> dAB_rec_aq(count_AB_rec, 1. - q);

    size_t ab_to_Ab = 0.5 * dab_rec_p(gen) + gaussianJitter_(gen);
    size_t ab_to_aB = 0.5 * dab_rec_q(gen) + gaussianJitter_(gen);
    size_t ab_stay_ab = count_ab_ - ab_to_Ab - ab_to_aB;

    size_t Ab_to_ab = 0.5 * dAb_rec_ap(gen) + gaussianJitter_(gen);
    size_t Ab_to_AB = 0.5 * dAb_rec_q(gen) + gaussianJitter_(gen);
    size_t Ab_stay_Ab = count_Ab_ - Ab_to_ab - Ab_to_AB;

    size_t aB_to_AB = 0.5 * daB_rec_p(gen) + gaussianJitter_(gen);
    size_t aB_to_ab = 0.5 * daB_rec_aq(gen) + gaussianJitter_(gen);
    size_t aB_stay_aB = count_aB_ - aB_to_AB - aB_to_ab;

    size_t AB_to_aB = 0.5 * dAB_rec_ap(gen) + gaussianJitter_(gen);
    size_t AB_to_Ab = 0.5 * dAB_rec_aq(gen) + gaussianJitter_(gen);
    size_t AB_stay_AB = count_AB_ - AB_to_aB - AB_to_Ab;

    count_ab_ = ab_stay_ab + Ab_to_ab + aB_to_ab;
    count_Ab_ = Ab_stay_Ab + ab_to_Ab + AB_to_Ab;
    count_aB_ = aB_stay_aB + ab_to_aB + AB_to_aB;
    count_AB_ = AB_stay_AB + Ab_to_AB + aB_to_AB;

    updateProps_();
  }

  void selectRandom_(std::mt19937& gen, double s)
  {
    assert(count_ab_ + count_Ab_ + count_aB_ + count_AB_ == n_);

    std::binomial_distribution<size_t> dAb(count_Ab_, std::abs(s));
    std::binomial_distribution<size_t> dAB(count_AB_, std::abs(s));

    size_t count_Ab_dead = dAb(gen);
    size_t count_AB_dead = dAB(gen);

    unsigned int Ab_replacement[4] = {0, 0, 0, 0};
    gsl_ran_multinomial(gen, n_, count_Ab_dead, { count_ab_, count_Ab_, count_aB_, count_AB_ }, Ab_replacement);

    unsigned int AB_replacement[4] = {0, 0, 0, 0};
    gsl_ran_multinomial(gen, n_, count_AB_dead, { count_ab_, count_Ab_, count_aB_, count_AB_ }, AB_replacement);

    count_ab_ = count_ab_ + Ab_replacement[0] + AB_replacement[0];
    count_Ab_ = count_Ab_ - count_Ab_dead + Ab_replacement[1] + AB_replacement[1];
    count_aB_ = count_aB_ + Ab_replacement[2] + AB_replacement[2];
    count_AB_ = count_AB_ - count_AB_dead + AB_replacement[3] + AB_replacement[3];

    delete [] Ab_replacement;
    delete [] AB_replacement;
  }

  void recombineDet_(double r)
  {
    prop_ab_ = prop_ab_ - r * fetchLD();
    prop_Ab_ = prop_Ab_ + r * fetchLD();
    prop_aB_ = prop_aB_ + r * fetchLD();
    prop_AB_ = prop_AB_ - r * fetchLD();

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

