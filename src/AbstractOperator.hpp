/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 30/05/2024
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
#include <algorithm>
#include <numeric>

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/unsupported/Eigen/MatrixFunctions>

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
  // flexible vector: one matrix per population (Drift, Mutation, Recombination and Selection) or pair thereof (Migration, Admixture)
  // the overal strategy is that matrices_ are built with coefficients only, and assigned indices that depend on the number of populations
  // they are then multiplied by parameters (1/N_i for Drift, m_ij for Migration etc) and finally added into transition_
  // this way the matrices_ need not be rebuilt during optimization when parameters change (see updateMatrices_() inside each derived class)
  std::vector<Eigen::SparseMatrix<long double>> matrices_; // "delta" matrix(ces)
  Eigen::SparseMatrix<long double> identity_; // helper matrix to convert from "delta" to "transition" matrix
  Eigen::SparseMatrix<long double> transition_; // "transition" matrix
  bpp::ParameterList prevParams_; // params in immediately previous iteration of optimization (for fast matrix updates)
  std::vector<size_t> popIndices_;

public:
  AbstractOperator():
  bpp::AbstractParameterAliasable(""),
  matrices_(0),
  identity_(),
  transition_(),
  prevParams_(),
  popIndices_(0)
  { }

  AbstractOperator(const std::vector<size_t>& popIndices):
  bpp::AbstractParameterAliasable(""),
  matrices_(0),
  identity_(),
  transition_(),
  prevParams_(),
  popIndices_(popIndices)
  { }

public:
  virtual ~AbstractOperator()
  {
    std::vector<std::string> paramNames(0);
    paramNames.reserve(getParameters().size());

    for(size_t i = 0; i < getParameters().size(); ++i)
      paramNames.emplace_back(getParameters()[i].getName());

    deleteParameters_(paramNames);
  }

  void setParameters(const bpp::ParameterList& params)
  {
    bpp::AbstractParameterAliasable::setParametersValues(params);
  }

  void fireParameterChanged(const bpp::ParameterList& params)
  {
    if(matchParametersValues(params))
      updateMatrices_();
  }

  const std::vector<Eigen::SparseMatrix<long double>>& getMatrices()
  {
    return matrices_; // delta matrices
  }

  const Eigen::SparseMatrix<long double>& getMatrix(size_t index)
  {
    return matrices_[index]; // delta matrix; population index for Drift, population-pair index for Migration etc
  }

  const Eigen::SparseMatrix<long double>& getTransitionMatrix()
  {
    return transition_;
  }

  const Eigen::SparseMatrix<long double>& getIdentity()
  {
    return identity_;
  }

  virtual void printDeltaLDMat(const std::string& fileName);

  virtual void printTransitionLDMat(const std::string& fileName);

protected:
  // this method sets up so-called "delta" matrices which govern the *change* in Y due to the operator
  virtual void setUpMatrices_(const SumStatsLibrary& sslib) = 0;  // called only once in order to set the coefficients

  virtual void updateMatrices_() = 0; // scales coefficients of "delta" matrices by (new) parameters during optimization

  // adds together the different matrices that make up an operator (one per population for Drift; population-pair for Migration, etc)
  virtual void assembleTransitionMatrix_();

  void setIdentity_(size_t numStats);

};

#endif
