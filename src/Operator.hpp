/*
 * Authors: Gustavo V. Barroso
 * Created: 29/07/2022
 * Last modified: 08/09/2022
 *
 */


#ifndef _OPERATOR_H_
#define _OPERATOR_H_

#include <vector>
#include <string>
#include <iostream>
#include <cmath>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <Bpp/Numeric/Function/Functions.h>
#include <Bpp/Numeric/AbstractParameterAliasable.h>
#include <Bpp/Numeric/Constraints.h>
#include <Bpp/Numeric/Parameter.h>
#include <Bpp/Numeric/ParameterList.h>
#include <Bpp/Text/TextTools.h>

#include "SumStatsLibrary.hpp"
#include "EigenDecomposition.hpp"

class Operator:
  public bpp::AbstractParameterAliasable
{

protected:
  // flexible vector: one matrix per population (Drift) or pair thereof (Migration), single matrices for Mutation and Recombination (?) etc
  // the overal strategy is that matrices_ are built with coefficients only, and assigned indices that depend on the number of populations
  // they are then multiplied by parameters (1/N_i for Drift, m_ij for Migration etc) and finally added into combinedMatrix_
  // this way the matrices_ need not be rebuilt during optimization when parameters change (see update_() inside each derived class)
  std::vector<Eigen::MatrixXd> matrices_; // NOTE do we want to keep these?
  std::vector<EigenDecomposition> eigenDec_;
  bpp::ParameterList prevParams_; // params in immediately previous iteration of optimization
  size_t exponent_; // for matrix exponentiation

public:
  Operator():
  bpp::AbstractParameterAliasable(""),
  matrices_(0),
  eigenDec_(0),
  prevParams_(),
  exponent_(1)
  { }

public:
  ~Operator()
  {
    std::cout << "Destruction of Operator with parameters:\n";
    getParameters().printParameters(std::cout);
  }

  Operator* clone() const
  {
    return new Operator(*this);
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

  const std::vector<Eigen::MatrixXd>& getMatrices()
  {
    return matrices_;
  }

  const std::vector<EigenDecomposition>& getEigenDecompositions()
  {
    return eigenDec_;
  }

  size_t getExponent()
  {
    return exponent_;
  }

  void setExponent(size_t exponent)
  {
    exponent_ = exponent;
  }

  const Eigen::MatrixXd& getMatrix(size_t index)
  {
    // population index for Drift; population-pair index for Migration etc
    return matrices_[index];
  }

  // adds together the different matrices that make up an operator (one per population for Drift; population-pair for Migration, etc)
  Eigen::MatrixXd fetchCombinedMatrix()
  {
    // NOTE that matrix exponentiation is built-in the Eigen Decomposition by construction of setUpMatrices_() and updateMatrices_()
    Eigen::MatrixXd mat = eigenDec_[0].matrixInverse() * eigenDec_[0].lambdaReal() * eigenDec_[0].matrix();

    if(eigenDec_.size() > 1)
    {
      for(size_t i = 1; i < eigenDec_.size(); ++i)
        mat += eigenDec_[i].matrixInverse() * eigenDec_[i].lambdaReal() * eigenDec_[i].matrix();
    }

    // adding Identity to convert from "delta" to "transition" matrix WARNING
    for(long int j = 0; j < mat.cols(); ++j)
      mat(j, j) += 1.;

    return mat;
  }

protected:
  virtual void setUpMatrices_(const SumStatsLibrary& sslib);  // called only once in order to set the coefficients

  virtual void updateMatrices_(); // scales coefficients by (new) parameters during optimization

};

#endif
