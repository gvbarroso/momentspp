/*
 * Authors: Gustavo V. Barroso
 * Created:29/07/2022
 * Last modified: 23/08/2022
 *
 */


#ifndef _OPERATOR_H_
#define _OPERATOR_H_

#include <vector>
#include <string>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Dense>

#include <Bpp/Numeric/Function/Functions.h>
#include <Bpp/Numeric/AbstractParameterAliasable.h>
#include <Bpp/Numeric/Constraints.h>
#include <Bpp/Numeric/Parameter.h>
#include <Bpp/Numeric/ParameterList.h>

#include "SumStatsLibrary.hpp"

class Operator:
  public bpp::AbstractParameterAliasable
{

private:
  // flexible vector: one matrix per population (Drift) or pair thereof (Migration), single matrices for Mutation and Recombination (?) etc
  // the overal strategy is that matrices_ are built with coefficients only, and assigned indices that depend on the number of populations
  // they are then multiplied by parameters (1/N_i for Drift, m_ij for Migration etc) and finally combined into combinedMatrix_
  // this way the matrices_ need not be rebuilt during optimization when parameters change (see update_() in each derived class)
  std::vector<Eigen::SparseMatrix<double>> matrices_;
  Eigen::SparseMatrix<double, Dynamic, Dynamic> combinedMatrix_;

  bpp::ParameterList prevParams_; // parameters values in immediately previous iteration of optimization

public:
  Operator():
  AbstractParameterAliasable(""),
  matrices_(0),
  combinedMatrix_(),
  prevParams_(0)
  { }

public:
  Operator* clone() const
  {
    return new Operator(*this);
  }

  void setParameters(const bpp::ParameterList& params)
  {
    AbstractParameterAliasable::setParametersValues(params);
  }

  void Operator::fireParameterChanged(const bpp::ParameterList& params)
  {
    matchParametersValues(params);
    updateMatrices();
  }

  const Eigen::SparseMatrix<int, Dynamic, Dynamic>& getMatrices()
  {
    return matrices_;
  }

  const Eigen::SparseMatrix<int, Dynamic, Dynamic>& getMatrix(size_t index)
  {
    // population index for Drift; population-pair index for Migration etc
    return matrices_[index];
  }

  const Eigen::SparseMatrix<int, Dynamic, Dynamic>& getCombinedMatrix()
  {
    return combinedMatrix_;
  }

   // scales matrix(ces) coefficients by (new) parameters values
  virtual void updateMatrices();

private:
  // setUpMatrices_ is meant to be called only once for each operator, in order to set the matrix(ces) coefficients
  virtual void setUpMatrices_(const SumStatsLibrary& sslib);

  void compressSparseMatrices_() // from Eigen perspective, NOT w.r.t moments with the same expectation
  {
    if(matrices_.size() == 0)
      throw bpp::Exception("Operator::tried to compress non-existing matrices!");

    else
      for(size_t i = 0; i < matrices_.size(); ++i)
        matrices_[i].makeCompressed();
  }

  void combineMatrices_()
  {
    combinedMatrix_ = matrices_[0];

    if(matrices_.size() > 1)
    {
      for(size_t i = 1; i <= matrices_.size(); ++i)
        combinedMatrix_ += matrices_[i];
    }
  }

};

#endif
