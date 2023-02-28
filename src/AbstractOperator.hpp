/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 21/10/2022
 *
 */


#ifndef _OPERATOR_H_
#define _OPERATOR_H_

#include <ios>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <memory>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <Bpp/Numeric/Function/Functions.h>
#include <Bpp/Numeric/AbstractParameterAliasable.h>
#include <Bpp/Numeric/Constraints.h>
#include <Bpp/Numeric/Parameter.h>
#include <Bpp/Numeric/ParameterList.h>
#include <Bpp/Text/TextTools.h>

#include "SumStatsLibrary.hpp"
#include "Log.hpp"

class AbstractOperator: public bpp::AbstractParameterAliasable
{

protected:
  // flexible vector: one matrix per population (Drift) or pair thereof (Migration), single matrices for Mutation and Recombination (?) etc
  // the overal strategy is that matrices_ are built with coefficients only, and assigned indices that depend on the number of populations
  // they are then multiplied by parameters (1/N_i for Drift, m_ij for Migration etc) and finally added into combinedMatrix_
  // this way the matrices_ need not be rebuilt during optimization when parameters change (see update_() inside each derived class)
  std::vector<Eigen::SparseMatrix<double>> matrices_;
  Eigen::SparseMatrix<double> identity_; // helper matrix to convert from "delta" to "transition" matrix

  // NOTE have combined matrix as a member to avoid assembling it when *this params don;t change in optimization?

  bpp::ParameterList prevParams_; // params in immediately previous iteration of optimization (for fast matrix updates)

public:
  AbstractOperator(size_t numStats):
  bpp::AbstractParameterAliasable(""),
  matrices_(0),
  identity_(),
  prevParams_()
  {
    setIdentity_(numStats);
  }

public:
  virtual ~AbstractOperator()
  { }

  void setParameters(const bpp::ParameterList& params)
  {
    bpp::AbstractParameterAliasable::setParametersValues(params);
  }

  void fireParameterChanged(const bpp::ParameterList& params)
  {
    if(matchParametersValues(params))
      updateMatrices_();
  }

  const std::vector<Eigen::SparseMatrix<double>>& getMatrices()
  {
    return matrices_;
  }

  const Eigen::SparseMatrix<double>& getMatrix(size_t index)
  {
    return matrices_[index]; // population index for Drift; population-pair index for Migration etc
  }

  // adds together the different matrices that make up an operator (one per population for Drift; population-pair for Migration, etc)
  virtual Eigen::SparseMatrix<double> fetchCombinedMatrix()
  {
    Eigen::SparseMatrix<double> mat = matrices_[0];

    if(matrices_.size() > 1)
    {
      for(size_t i = 1; i < matrices_.size(); ++i)
        mat += matrices_[i];
    }

    return mat += identity_; // adds Identity to convert from "delta" to "transition" matrix (fixes row of Dummy moment)
  }

public:
  // this method sets up so-called "delta" matrices which govern the *change* in Y due to the operator
  virtual void setUpMatrices_(const SumStatsLibrary& sslib) = 0;  // called only once in order to set the coefficients

  virtual void updateMatrices_() = 0; // scales coefficients of "delta" matrices by (new) parameters during optimization

  void setIdentity_(size_t numStats)
  {
    Eigen::SparseMatrix<double> id(numStats, numStats);

    std::vector<Eigen::Triplet<double>> md(0);
    md.reserve(numStats);

    for(size_t i = 0; i < numStats; ++i)
      md.emplace_back(Eigen::Triplet<double>(i, i, 1.));

    id.setFromTriplets(std::begin(md), std::end(md));
    id.makeCompressed();

    identity_ = id;
  }

};

#endif
