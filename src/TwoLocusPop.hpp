/*
 * Authors: Gustavo V. Barroso
 * Created: 24/10/2023
 * Last modified: 28/10/2024
 */


#ifndef _TWO_LOCUS_POP_H_
#define _TWO_LOCUS_POP_H_

#include <algorithm>
#include <utility>
#include <vector>
#include <iostream>
#include <array>
#include <numeric>

#include <boost/algorithm/string.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/zlib.hpp>

class TwoLocusPop
{

private:
  size_t id_;
  size_t l_; // length of "segment"
  unsigned int ne_; // pop size
  std::vector<std::array<double, 4>> x_; // in order: fAB, fAb, faB, fab
  std::vector<double> xl_; // allele freqs. at (independent) left locus
  std::vector<double> xr_; // allele freqs. at (independent) right locus

  std::vector<double> hl_;
  std::vector<double> hr_;
  std::vector<double> d_;
  std::vector<double> dr_; // D(1-2q)
  std::vector<double> dl_; // D(1-2p)
  std::vector<double> dz_; // D(1-2p)(1-2q)
  std::vector<double> dsqr_;
  std::vector<double> pi2_;

public:
  TwoLocusPop():
  id_(0),
  l_(1),
  ne_(1000),
  x_(0),
  xl_(0),
  xr_(0),
  hl_(0),
  hr_(0),
  d_(0),
  dr_(0),
  dl_(0),
  dz_(0),
  dsqr_(0),
  pi2_(0)
  { }

  TwoLocusPop(size_t id, size_t l, unsigned int ne):
  id_(id),
  l_(l),
  ne_(ne),
  x_(0),
  xl_(0),
  xr_(0),
  hl_(0),
  hr_(0),
  d_(0),
  dr_(0),
  dl_(0),
  dz_(0),
  dsqr_(0),
  pi2_(0)
  {
    x_.reserve(10 * l);
    xl_.reserve(10 * l);
    xr_.reserve(10 * l);
    d_.reserve(10 * l);
    dr_.reserve(10 * l);
    dl_.reserve(10 * l);
    dz_.reserve(10 * l);
    dsqr_.reserve(10 * l);
    hl_.reserve(10 * l);
    hr_.reserve(10 * l);
    pi2_.reserve(10 * l);
  }

public:
  void printXl(std::ostream& stream)
  {
    stream << "Xl:\n\t";

    for(auto& v : xl_)
      stream << v << ",";

    stream << "\n";
  }

  void printXr(std::ostream& stream)
  {
    stream << "Xr:\n\t";

    for(auto& v : xr_)
      stream << v << ",";

    stream << "\n";
  }

  void printX(std::ostream& stream)
  {
    for(size_t i = 0; i < x_.size(); ++i)
    {
      stream << x_[i][0] << "\t";
      stream << x_[i][1] << "\t";
      stream << x_[i][2] << "\t";
      stream << x_[i][3];
    }

    stream << "\n";
  }

  void printAttributes(std::ostream& stream)
  {
    stream << "Pop. id: " << id_ << "\n";

    printXl(stream);
    printXr(stream);
    printX(stream);
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

  void compute_Hls()
  {
    hl_.resize(xl_.size());

    for(size_t i = 0; i < xl_.size(); ++i)
      hl_[i] = xl_[i] * (1. - xl_[i]);
  }

  void compute_Hrs()
  {
    hr_.resize(xr_.size());

    for(size_t i = 0; i < xr_.size(); ++i)
      hr_[i] = xr_[i] * (1. - xr_[i]);
  }

  void compute_Ds()
  {
    d_.resize(x_.size());

    for(size_t i = 0; i < x_.size(); ++i)
      d_[i] = x_[i][0] * x_[i][3] - x_[i][1] * x_[i][2];
  }

  void compute_Drs()
  {
    dr_.resize(x_.size());

    for(size_t i = 0; i < x_.size(); ++i)
      dr_[i] = (x_[i][0] * x_[i][3] - x_[i][1] * x_[i][2]) * (1 - 2 * (x_[i][0] + x_[i][2]));
  }

  void compute_Dls()
  {
    dl_.resize(x_.size());

    for(size_t i = 0; i < x_.size(); ++i)
      dl_[i] = (x_[i][0] * x_[i][3] - x_[i][1] * x_[i][2]) * (1 - 2 * (x_[i][0] + x_[i][1]));
  }

  void compute_Dzs()
  {
    dz_.resize(x_.size());

    for(size_t i = 0; i < x_.size(); ++i)
    {
      double d = x_[i][0] * x_[i][3] - x_[i][1] * x_[i][2];
      double p = x_[i][0] + x_[i][1];
      double q = x_[i][0] + x_[i][2];

      dz_[i] = d * (1. - 2. * p) * (1. - 2. * q);
    }
  }

  void compute_Dsqrs()
  {
    dsqr_.resize(x_.size());

    for(size_t i = 0; i < x_.size(); ++i)
      dsqr_[i] = (x_[i][0] * x_[i][3] - x_[i][1] * x_[i][2]) * (x_[i][0] * x_[i][3] - x_[i][1] * x_[i][2]);
  }

  void compute_Pi2s()
  {
    pi2_.resize(x_.size());

    for(size_t i = 0; i < x_.size(); ++i)
    {
      double p = x_[i][0] + x_[i][1];
      double q = x_[i][0] + x_[i][2];

      pi2_[i] = p * (1. - p) * q * (1. - q);
    }
  }

  double getSumHl()
  {
    std::sort(std::begin(hl_), std::end(hl_));
    return std::accumulate(std::begin(hl_), std::end(hl_), 0.);
  }

  double getSumHr()
  {
    std::sort(std::begin(hr_), std::end(hr_));
    return std::accumulate(std::begin(hr_), std::end(hr_), 0.);
  }

  double getSumD()
  {
    std::sort(std::begin(d_), std::end(d_));
    return std::accumulate(std::begin(d_), std::end(d_), 0.);
  }

  double getSumDr()
  {
    std::sort(std::begin(dr_), std::end(dr_));
    return std::accumulate(std::begin(dr_), std::end(dr_), 0.);
  }

  double getSumDl()
  {
    std::sort(std::begin(dl_), std::end(dl_));
    return std::accumulate(std::begin(dl_), std::end(dl_), 0.);
  }

  double getSumDz()
  {
    std::sort(std::begin(dz_), std::end(dz_));
    return std::accumulate(std::begin(dz_), std::end(dz_), 0.);
  }

  double getSumDsqr()
  {
    std::sort(std::begin(dsqr_), std::end(dsqr_));
    return std::accumulate(std::begin(dsqr_), std::end(dsqr_), 0.);
  }

  double getSumPi2()
  {
    std::sort(std::begin(pi2_), std::end(pi2_));
    return std::accumulate(std::begin(pi2_), std::end(pi2_), 0.);
  }

  void computeStats()
  {
    compute_Hls();
    compute_Hrs();
    compute_Ds();
    compute_Drs();
    compute_Dls();
    compute_Dzs();
    compute_Dsqrs();
    compute_Pi2s();
  }

  void tabulate_Ds(boost::iostreams::filtering_ostream& file_d)
  {
    for(size_t i = 0; i < x_.size(); ++i)
    {
      double d = x_[i][0] * x_[i][3] - x_[i][1] * x_[i][2];
      double p = x_[i][0] + x_[i][1];
      double q = x_[i][0] + x_[i][2];

      file_d << d << "\t" << p << "\t" << q << "\n";
    }
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
      double p = x_[i][0] + x_[i][1];

      x_[i][0] *= (1. + s * (1. - p));
      x_[i][1] *= (1. + s * (1. - p));
      x_[i][2] *= (1. - s * p);
      x_[i][3] *= (1. - s * p);
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

    // new (independent) mutations at left locus
    unsigned int num_left_mut = gsl_ran_poisson(gen, 2. * ne_ * u * l_);
    for(size_t i = 0; i < num_left_mut; ++i)
      xl_.emplace_back(1. / (2. * ne_));

    // new (independent) mutations at right locus
    unsigned int num_right_mut = gsl_ran_poisson(gen, 2. * ne_ * u * l_);
    for(size_t i = 0; i < num_right_mut; ++i)
      xr_.emplace_back(1. / (2. * ne_));
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

  void cleanup(std::ostream& log)
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
      {
        it = x_.erase(it);
        log << "removed\ttwo-locus\tsystem\tnum. pairs = " << x_.size() << "\n";
      }
      else
        ++it;
    }
  }

private:
  void migrate_(TwoLocusPop& other, double m);

};

#endif

