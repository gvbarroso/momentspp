/*
 * Authors: Gustavo V. Barroso
 * Created: 10/08/2022
 * Last modified: 26/09/2022
 *
 */


#ifndef _MUTATION_H_
#define _MUTATION_H_

#include "AbstractOperator.hpp"

class Mutation:
  public AbstractOperator
{

protected:
  Eigen::SparseMatrix<double> oneLocusPi_;

public:
  Mutation(std::shared_ptr<bpp::IntervalConstraint> ic, const SumStatsLibrary& ssl):
  AbstractOperator(),
  oneLocusPi_()
  {
    addParameter_(new bpp::Parameter("u_0", 1e-8, ic));

    prevParams_.addParameters(getParameters()); // inits list of "previous" parameters
    setUpMatrices_(ssl);
  }

  virtual Mutation* clone() const override
  {
    return new Mutation(*this);
  }

  const Eigen::SparseMatrix<double>& getOneLocusPi() const
  {
    return oneLocusPi_;
  }

  void setUpMatrices_(const SumStatsLibrary& ssl);

  void updateMatrices_();

};

#endif
