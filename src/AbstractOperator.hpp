/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 12/09/2025
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

#include <omp.h>

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/unsupported/Eigen/MatrixFunctions>
#include <eigen3/unsupported/Eigen/MPRealSupport> // for arbitrary-precision arithmetic

#include <Bpp/Numeric/Function/Functions.h>
#include <Bpp/Numeric/AbstractParameterAliasable.h>
#include <Bpp/Numeric/Constraints.h>
#include <Bpp/Numeric/Parameter.h>
#include <Bpp/Numeric/ParameterList.h>
#include <Bpp/Text/TextTools.h>

#include "MatrixEngine.hpp"
#include "SumStatsLibrary.hpp"
#include "Log.hpp"

class AbstractOperator: public bpp::AbstractParameterAliasable
{

protected:
  // flexible vector: one matrix per population (Drift, Mutation, Recombination and Selection) or pair thereof (Migration, Admixture)
  // the overal strategy is that matrices_ are built with coefficients only, and assigned indices that depend on the number of populations
  // they are then multiplied by parameters (1/2N_i for Drift, m_ij for Migration etc) and finally added into transition_
  // this way the matrices_ need not be rebuilt during optimization when parameters change (see updateMatrices_() inside each derived class)
  std::vector<std::unique_ptr<MatrixEngine>> matrices_; // "delta" matrix(ces)
  std::unique_ptr<MatrixEngine> transition_;  // "transition" matrix
  MatrixVariant identityMatrix_;
  bool identityInitialized_;
  bpp::ParameterList prevParams_; // params in immediately previous iteration of optimization (for fast matrix updates)
  std::vector<size_t> popIndices_;

public:
  AbstractOperator():
  bpp::AbstractParameterAliasable(""),
  matrices_(0),
  transition_(nullptr),
  identityMatrix_(),
  identityInitialized_(false),
  prevParams_(),
  popIndices_(0)
  { }

  AbstractOperator(const std::vector<size_t>& popIndices):
  bpp::AbstractParameterAliasable(""),
  matrices_(0),
  transition_(nullptr),
  identityMatrix_(),
  identityInitialized_(false),
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

  const std::vector<std::unique_ptr<MatrixEngine>>& getMatrices()
  {
    return matrices_;
  }

  const std::unique_ptr<MatrixEngine>& getMatrix(size_t index)
  {
    return matrices_[index];
  }

  const std::unique_ptr<MatrixEngine>& getTransitionMatrix()
  {
    return transition_;
  }

  virtual void printDeltaLDMat(const std::string& fileName);

protected:
  // this method sets up so-called "delta" matrices which govern the *change* in Y due to the operator
  virtual void setUpMatrices_(const SumStatsLibrary& sslib, bool highPrecision) = 0;  // called only once in order to set the coefficients

  virtual void updateMatrices_() = 0; // scales coefficients of "delta" matrices by (new) parameters during optimization

  // adds together the different matrices that make up an operator (one per population for Drift; population-pair for Migration, etc)
  virtual void assembleTransitionMatrix_();

};

#endif
