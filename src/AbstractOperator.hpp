/*
 * Authors: Gustavo V. Barroso
 * Created: 09/08/2022
 * Last modified: 24/09/2025
 *
 */

#ifndef _OPERATOR_H_
#define _OPERATOR_H_

#include "VariantUtils.hpp"
#include "MatrixEngine.hpp"
#include "SumStatsLibrary.hpp"
#include "Log.hpp"

#include <Bpp/Numeric/AbstractParameterAliasable.h>
#include <Bpp/Numeric/ParameterList.h>

#include <vector>
#include <memory>
#include <string>


class AbstractOperator : public bpp::AbstractParameterAliasable
{
protected:
  // flexible vector: one matrix per population (Drift, Mutation, Recombination and Selection) or pair thereof (Migration, Admixture)
  // the overal strategy is that matrices_ are built with coefficients only, and assigned indices that depend on the number of populations
  // they are then multiplied by parameters (1/N_i for Drift, m_ij for Migration etc) and finally added into transition_
  // this way the matrices_ need not be rebuilt during optimization when parameters change (see updateMatrices_() inside each derived class)
  std::vector<std::unique_ptr<MatrixEngine>> matrices_; // "delta" matrix(ces)
  std::unique_ptr<MatrixEngine> transition_; // "transition" matrix in Sparse format (sum of all `matrices_` after scaling by current params_
  bpp::ParameterList prevParams_; // params in immediately previous iteration of optimization (for fast matrix updates)
  std::vector<size_t> popIndices_; // which populations *this operator acts on

public:
  AbstractOperator() noexcept:
  bpp::AbstractParameterAliasable(""),
  matrices_(),
  transition_(nullptr),
  prevParams_(),
  popIndices_()
  { }

  explicit AbstractOperator(const std::vector<size_t>& popIndices, bool highPrecision) noexcept:
  bpp::AbstractParameterAliasable(""),
  matrices_(),
  transition_(nullptr),
  prevParams_(),
  popIndices_(popIndices)
  {
    initializeScalarType_(highPrecision);
  }

  // Deep-copy: clones each MatrixEngine
  AbstractOperator(const AbstractOperator& other):
  bpp::AbstractParameterAliasable(""),
  prevParams_(other.prevParams_),
  popIndices_(other.popIndices_)
  {
    matrices_.clear();
    matrices_.reserve(other.matrices_.size());
    for(auto const& matPtr : other.matrices_)
      matrices_.emplace_back(matPtr ? matPtr->clone() : nullptr);

    transition_ = other.transition_ ? other.transition_->clone() : nullptr;
  }

  AbstractOperator(AbstractOperator&&) noexcept = default;

  AbstractOperator& operator=(const AbstractOperator& other)
  {
    if(this != &other)
    {
      matrices_.clear();
      matrices_.reserve(other.matrices_.size());

      for(auto const& matPtr : other.matrices_)
        matrices_.emplace_back(matPtr ? matPtr->clone() : nullptr);

      transition_ = other.transition_ ? other.transition_->clone() : nullptr;
      prevParams_ = other.prevParams_;
      popIndices_ = other.popIndices_;
    }

    return *this;
  }

  // named swap for strong exception safety
  void swap(AbstractOperator& other) noexcept
  {
    std::swap(matrices_, other.matrices_);
    std::swap(transition_, other.transition_);
    std::swap(prevParams_, other.prevParams_);
    std::swap(popIndices_, other.popIndices_);
    // base class swap is not needed
  }

  friend void swap(AbstractOperator& a, AbstractOperator& b) noexcept
  {
    a.swap(b);
  }

  // polymorphic clone
  AbstractOperator* clone() const override = 0;

  std::unique_ptr<AbstractOperator> cloneOperator() const
  {
    return std::unique_ptr<AbstractOperator>(clone());
  }

  virtual ~AbstractOperator()
  {
    std::vector<std::string> paramNames(0);
    paramNames.reserve(getParameters().size());

    for(size_t i = 0; i < getParameters().size(); ++i)
      paramNames.emplace_back(getParameters()[i].getName());

    deleteParameters_(paramNames);
  }

  // Replace entire parameter set; triggers update if any changed
  void setParameters(const bpp::ParameterList& params)
  {
    bpp::AbstractParameterAliasable::setParametersValues(params);
  }

  void fireParameterChanged(const bpp::ParameterList& params) override
  {
    if(matchParametersValues(params))
      updateMatrices_();
  }

  const std::vector<size_t>& getPopIndices() const noexcept
  {
    return popIndices_;
  }

  const std::vector<std::unique_ptr<MatrixEngine>>& getMatrices() const noexcept
  {
    return matrices_;
  }

  const MatrixEngine& getMatrix(size_t idx) const
  {
    return *matrices_.at(idx);
  }

  // extracts the raw MatrixVariant — empty if no transition set
  MatrixEngine::MatrixVariant getTransitionMatrixVariant() const
  {
    if(!transition_)
      return MatrixEngine::MatrixVariant{};

    return transition_->getMatrixVariant();
  }

  // extracts the Eigen‐view variant (dense/sparse) — empty if no transition
  MatrixEngine::MatrixEigenVariant getTransitionMatrixVariantEigen() const
  {
    if(!transition_)
      return MatrixEngine::MatrixEigenVariant{};

    return transition_->toEigenMatrixVariant();
  }

  // mostly for debugging
  virtual void printDeltaLDMat(const std::string& fileName);

  // scales a matrix of coefficients (eg by parameter value)
  void scaleMatrix(double factor)
  {
    if(!transition_)
      throw bpp::Exception("scaleMatrix: transition matrix is not set");

    transition_->scaleMatrix(factor);
  }

protected:
  // so that setUpMatrices_() knows which Scalar type to use
  void initializeScalarType_(bool highPrecision)
  {
    if(highPrecision)
      transition_ = MatrixEngine::createEmpty<mpfr::mpreal>(1);

    else
      transition_ = MatrixEngine::createEmpty<double>(1);
  }

  virtual void setUpMatrices_(const SumStatsLibrary& sslib) = 0;

  virtual void updateMatrices_() = 0;

  virtual void assembleTransitionMatrix_()
  {
    if(matrices_.empty())
    {
      transition_.reset();
      return;
    }

    // clones first delta as the base
    transition_ = matrices_[0]
      ? matrices_[0]->clone()
      : nullptr;

    // adds all others via your variant‐safe helper
    for(size_t i = 1; i < matrices_.size(); ++i)
    {
      if(!matrices_[i])
        continue;

      transition_->addToMatrix(matrices_[i]->getMatrixVariant());
    }
  }

};

#endif
