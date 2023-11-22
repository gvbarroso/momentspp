/*
 * Authors: Gustavo V. Barroso
 * Created: 24/10/2023
 * Last modified: 22/11/2023
 */


#ifndef _TWO_LOCUS_POP_H_
#define _TWO_LOCUS_POP_H_

#include <algorithm>
#include <utility>
#include <vector>
#include <iostream>
#include <array>
#include <numeric>

class TwoLocusPop
{

private:
  size_t id_;
  size_t l_; // length of segment
  unsigned int ne_; // pop size
  std::vector<std::array<double, 4>> x_; // in order: AB, Ab, aB, ab
  std::vector<double> xl_; // allele freqs. at (independent) left locus
  std::vector<double> xr_; // allele freqs. at (independent) right locus

public:
  TwoLocusPop():
  id_(0),
  l_(1),
  ne_(1000),
  x_(0),
  xl_(0),
  xr_(0)
  { }

  TwoLocusPop(size_t id, size_t l, unsigned int ne):
  id_(id),
  l_(l),
  ne_(ne),
  x_(0),
  xl_(0),
  xr_(0)
  {
    x_.reserve(10 * l);
    xl_.reserve(10 * l);
    xr_.reserve(10 * l);
  }

public:
  void printXl()
  {
    std::cout << "Xl:\n\t";

    for(auto& v : xl_)
      std::cout << v << ",";

    std::cout << "\n";
  }

  void printXr()
  {
    std::cout << "Xr:\n\t";

    for(auto& v : xr_)
      std::cout << v << ",";

    std::cout << "\n";
  }

  void printX()
  {
    for(size_t i = 0; i < x_.size(); ++i)
    {
      std::cout << "\tfAB = " << x_[i][0] << ", ";
      std::cout << "fAb = " << x_[i][1] << ", ";
      std::cout << "faB = " << x_[i][2] << ", ";
      std::cout << "fab = " << x_[i][3] << "\n";
    }

    std::cout << "\n";
  }

  void printAttributes(std::ostream& stream)
  {
    stream << "Pop. id: " << id_ << "\n";

    printXl();
    printXr();
    printX();
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

  void setPopSize(unsigned int ne)
  {
    ne_ = ne;
  }

  const std::vector<std::array<double, 4>>& getX() const
  {
    return x_;
  }

  std::vector<std::array<double, 4>>& getX()
  {
    return x_;
  }

  const std::vector<double>& getXl() const
  {
    return xl_;
  }

  std::vector<double>& getXl()
  {
    return xl_;
  }

  const std::vector<double>& getXr() const
  {
    return xr_;
  }

  std::vector<double>& getXr()
  {
    return xr_;
  }

  std::vector<double> getDs()
  {
    std::vector<double> ret(0);
    ret.reserve(x_.size());

    for(size_t i = 0; i < x_.size(); ++i)
      ret.emplace_back(x_[i][0] * x_[i][3] - x_[i][1] * x_[i][2]);

    return ret;
  }

  std::vector<double> getPs()
  {
    std::vector<double> ret(0);
    ret.reserve(x_.size());

    for(size_t i = 0; i < x_.size(); ++i)
      ret.emplace_back(x_[i][0] + x_[i][1]);

    return ret;
  }

  std::vector<double> getQs()
  {
    std::vector<double> ret(0);
    ret.reserve(x_.size());

    for(size_t i = 0; i < x_.size(); ++i)
      ret.emplace_back(x_[i][0] + x_[i][2]);

    return ret;
  }

  std::vector<double> getHls()
  {
    std::vector<double> ret(0);
    ret.reserve(xl_.size());

    for(size_t i = 0; i < xl_.size(); ++i)
      ret.emplace_back(xl_[i] * (1. - xl_[i]));

    return ret;
  }

  std::vector<double> getHrs()
  {
    std::vector<double> ret(0);
    ret.reserve(xr_.size());

    for(size_t i = 0; i < xr_.size(); ++i)
      ret.emplace_back(xr_[i] * (1. - xr_[i]));

    return ret;
  }

  std::vector<double> getDsqrs()
  {
    std::vector<double> ret(0);
    ret.reserve(x_.size());

    for(size_t i = 0; i < x_.size(); ++i)
      ret.emplace_back((x_[i][0] * x_[i][3] - x_[i][1] * x_[i][2]) * (x_[i][0] * x_[i][3] - x_[i][1] * x_[i][2]));

    return ret;
  }

  std::vector<double> getDzs()
  {
    std::vector<double> ret(0);
    ret.reserve(x_.size());

    for(size_t i = 0; i < x_.size(); ++i)
    {
      double d = x_[i][0] * x_[i][3] - x_[i][1] * x_[i][2];
      double p = x_[i][0] + x_[i][1];
      double q = x_[i][0] + x_[i][2];

      ret.emplace_back(d * (1. - 2. * p) * (1. - 2. * q));
    }

    return ret;
  }

  std::vector<double> getPi2s()
  {
    std::vector<double> ret(0);
    ret.reserve(x_.size());

    for(size_t i = 0; i < x_.size(); ++i)
    {
      double p = x_[i][0] + x_[i][1];
      double q = x_[i][0] + x_[i][2];

      ret.emplace_back(p * (1. - p) * q * (1. - q));
    }

    return ret;
  }

  double getSumHl()
  {
    std::vector<double> vals = getHls();
    std::sort(std::begin(vals), std::end(vals));
    return std::accumulate(std::begin(vals), std::end(vals), 0.);
  }

  double getSumHr()
  {
    std::vector<double> vals = getHrs();
    std::sort(std::begin(vals), std::end(vals));
    return std::accumulate(std::begin(vals), std::end(vals), 0.);
  }

  double getSumDsqr()
  {
    std::vector<double> vals = getDsqrs();
    std::sort(std::begin(vals), std::end(vals));
    return std::accumulate(std::begin(vals), std::end(vals), 0.);
  }

  double getSumDz()
  {
    std::vector<double> vals = getDzs();
    std::sort(std::begin(vals), std::end(vals));
    return std::accumulate(std::begin(vals), std::end(vals), 0.);
  }

  double getSumPi2()
  {
    std::vector<double> vals = getPi2s();
    std::sort(std::begin(vals), std::end(vals));
    return std::accumulate(std::begin(vals), std::end(vals), 0.);
  }

  void recombine(double r)
  {
    for(size_t i = 0; i < x_.size(); ++i)
    {
      double d = x_[i][0] * x_[i][3] - x_[i][1] * x_[i][2];

      x_[i][0] -= r * d;
      x_[i][1] += r * d;
      x_[i][2] += r * d;
      x_[i][3] -= r * d;
    }
  }

  void select(double s)
  {
    for(size_t i = 0; i < xl_.size(); ++i)
      xl_[i] = ((1 + s) * xl_[i]) / (xl_[i] * (1 + s) + 1 - xl_[i]);

    for(size_t i = 0; i < x_.size(); ++i)
    {
      //printX();
      //double before = x_[i][0] + x_[i][1] + x_[i][2] + x_[i][3];

      double p = x_[i][0] + x_[i][1];

      x_[i][0] *= (1. + s * (1. - p));
      x_[i][1] *= (1. + s * (1. - p));
      x_[i][2] *= (1. - s * p);
      x_[i][3] *= (1. - s * p);

      //double after = x_[i][0] + x_[i][1] + x_[i][2] + x_[i][3];

      //std::cout << std::setprecision(12) << "before = " << before << "; after = " << after << "\n";
      //printX();
    }
  }

  void drift(const gsl_rng* gen)
  {
    for(size_t i = 0; i < xl_.size(); ++i)
      xl_[i] = gsl_ran_binomial(gen, xl_[i], 2. * ne_) / (2. * ne_);

    for(size_t i = 0; i < xr_.size(); ++i)
      xr_[i] = gsl_ran_binomial(gen, xr_[i], 2. * ne_) / (2. * ne_);

    for(size_t i = 0; i < x_.size(); ++i)
    {
      unsigned int next[4];
      double probs[4] = { x_[i][0], x_[i][1], x_[i][2], x_[i][3] };

      gsl_ran_multinomial(gen, 4, 2 * ne_, probs, next);

      x_[i][0] = static_cast<double>(next[0]) / (2. * ne_);
      x_[i][1] = static_cast<double>(next[1]) / (2. * ne_);
      x_[i][2] = static_cast<double>(next[2]) / (2. * ne_);
      x_[i][3] = static_cast<double>(next[3]) / (2. * ne_);
    }
  }

  void mutate(double u, const gsl_rng* gen)
  {
    // new single mutations at left locus
    unsigned int num_left_mut = gsl_ran_poisson(gen, 2. * ne_ * u * l_);
    for(size_t i = 0; i < num_left_mut; ++i)
      xl_.emplace_back(1. / (2. * ne_));

    // new single mutations at right locus
    unsigned int num_right_mut = gsl_ran_poisson(gen, 2. * ne_ * u * l_);
    for(size_t i = 0; i < num_right_mut; ++i)
      xr_.emplace_back(1. / (2. * ne_));


    for(size_t i = 0; i < l_; ++i)
    {
      // mutate against segregating left loci in Xl
      for(size_t j = 0; j < xl_.size(); ++j)
      {
        if(gsl_rng_uniform(gen) < 2 * ne_ * u)
        {
          if(gsl_rng_uniform(gen) < xl_[j]) // falls on Ab background
          {
            std::array<double, 4> vals = {1. / (2. * ne_), xl_[j] - 1. / (2. * ne_), 0., 1. - xl_[j]};
            x_.emplace_back(vals);
          }

          else // falls on ab background
          {
            std::array<double, 4> vals = {0., xl_[j], 1. / (2. * ne_), 1. - xl_[j] - 1. / (2. * ne_)};
            x_.emplace_back(vals);
          }
        }
      }

      // mutate against segregating right loci in Xr
      for(size_t j = 0; j < xr_.size(); ++j)
      {
        if(gsl_rng_uniform(gen) < 2 * ne_ * u)
        {
          if(gsl_rng_uniform(gen) < xr_[j]) // falls on aB background
          {
            std::array<double, 4> vals = {1. / (2. * ne_), 0., xr_[j] - 1. / (2. * ne_), 1. - xr_[j]};
            x_.emplace_back(vals);
          }

          else // falls on ab background
          {
            std::array<double, 4> vals = {0., 1. / (2. * ne_), xr_[j], 1. - xr_[j] - 1. / (2. * ne_)};
            x_.emplace_back(vals);
          }
        }
      }
    }
  }

  void evolve(double u, double r, double s, const gsl_rng* gen)
  {
    recombine(r);
    select(s);
    drift(gen);
    mutate(u, gen);
  }

  void evolve(TwoLocusPop& other, double m, double u, double r, double s, const gsl_rng* gen)
  {
    migrate_(other, m);

    recombine(r);
    select(s);
    drift(gen);
    mutate(u, gen);
  }

  void cleanup()
  {
    // remove fixed or lost in xl_ and xr_
    for(auto it = std::begin(xl_); it != std::end(xl_);)
    {
      if((*it) == 0. || (*it) == 1.)
        it = xl_.erase(it);

      else
        ++it;
    }

    for(auto it = std::begin(xr_); it != std::end(xr_);)
    {
      if((*it) == 0. || (*it) == 1.)
        it = xr_.erase(it);

      else
        ++it;
    }

    // remove any pair with p or q at 0 or 1
    for(auto it = std::begin(x_); it != std::end(x_);)
    {
      double p = (*it)[0] + (*it)[1];
      double q = (*it)[0] + (*it)[2];

      if(p == 0. || p == 1. || q == 0. || q == 1.)
        it = x_.erase(it);

      else
        ++it;
    }
  }

private:
  void migrate_(TwoLocusPop& other, double m);

};

#endif

