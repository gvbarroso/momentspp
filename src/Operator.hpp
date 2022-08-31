/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 31/08/2022
 *
 */


#ifndef _OPERATOR_H_
#define _OPERATOR_H_

#include <vector>
#include <string>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>

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
  // they are then multiplied by parameters (1/N_i for Drift, m_ij for Migration etc) and finally added into combinedMatrix_
  // this way the matrices_ need not be rebuilt during optimization when parameters change (see update_() inside each derived class)
  std::vector<Eigen::SparseMatrix<double>> matrices_; // NOTE do we want to keep these?
  std::vector<Eigen::EigenSolver> solvers_;

  bpp::ParameterList prevParams_; // parameters values in immediately previous iteration of optimization

public:
  Operator():
  AbstractParameterAliasable(""),
  matrices_(0),
  solvers_(0),
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
    if(matchParametersValues(params))
      updateMatrices_();
  }

  const std::vector<Eigen::SparseMatrix<double>>& getMatrices()
  {
    return matrices_;
  }

  const std::vector<Eigen::EigenSolver<double>>& getEigenSolvers()
  {
    return solvers_;
  }

  const Eigen::SparseMatrix<double>& getMatrix(size_t index)
  {
    // population index for Drift; population-pair index for Migration etc
    return matrices_[index];
  }

  Eigen::SparseMatrix<double> fetchCombinedMatrix(size_t exponent)
  {
    // we want something like this: (TODO check if the EigenSolver returns const refs to these objects or computes them on the fly)
    Eigen::SparseMatrix<double> mat = solvers_[0].eigenvectors() * solvers_[0].eigenvalues() ^ exponent * solvers_[0].eigenvectors().inverse());

    if(matrices_.size() > 1)
      for(size_t i = 1; i < matrices_.size(); ++i)
        mat += solvers_[i].eigenvectors() * (solvers_[i].eigenvalues() ^ exponent).asDiagonal() * solvers_[i].eigenvectors().inverse())

    return mat;
  }

private:
  virtual void setUpMatrices_(const SumStatsLibrary& sslib);  // called only once in order to set the coefficients

  virtual void updateMatrices_(); // scales coefficients by (new) parameters during optimization

  void compressSparseMatrices_() // from Eigen3 perspective, NOT w.r.t moments with the same expectation
  {
    if(matrices_.size() == 0)
      throw bpp::Exception("Operator::tried to compress non-existing matrices!");

    else
      for(size_t i = 0; i < matrices_.size(); ++i)
        matrices_[i].makeCompressed();
  }

};

#endif
