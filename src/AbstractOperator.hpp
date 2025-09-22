#ifndef _OPERATOR_H_
#define _OPERATOR_H_

// ==== AbstractOperator.hpp ====
#pragma once

#include "VariantUtils.hpp"     // visitSameType + overloaded
#include "MatrixEngine.hpp"
#include "SumStatsLibrary.hpp"
#include "Log.hpp"

#include <Bpp/Numeric/AbstractParameterAliasable.h>
#include <Bpp/Numeric/ParameterList.h>

#include <vector>
#include <memory>
#include <string>

namespace bpp {
  class ParameterList;
  class Exception;
}

/// AbstractOperator drives one “step” of your pipeline by
/// building per-population (or per-pair) delta-matrices, then
/// assembling them into a single transition matrix.
class AbstractOperator : public bpp::AbstractParameterAliasable
{
protected:
  // One MatrixEngine per population or population-pair
  std::vector<std::unique_ptr<MatrixEngine>> matrices_;

  // Final sum of all `matrices_` after scaling by current params
  std::unique_ptr<MatrixEngine> transition_;

  // Track last parameters so updateMatrices_() knows when to rebuild
  bpp::ParameterList prevParams_;

  // Which populations this operator acts on
  std::vector<size_t> popIndices_;

public:
  // Default ctor: empty operator
  AbstractOperator() noexcept:
  bpp::AbstractParameterAliasable(""),
  matrices_(),
  transition_(nullptr),
  prevParams_(),
  popIndices_()
  { }

  // Main ctor: specify which populations to touch
  explicit AbstractOperator(const std::vector<size_t>& popIndices) noexcept:
  bpp::AbstractParameterAliasable(""),
  matrices_(),
  transition_(nullptr),
  prevParams_(),
  popIndices_(popIndices)
  { }

  // Deep-copy: clone each MatrixEngine
  AbstractOperator(const AbstractOperator& other):
  bpp::AbstractParameterAliasable(""),
  popIndices_(other.popIndices_),
  prevParams_(other.prevParams_)
  {
    matrices_.reserve(other.matrices_.size());
    for(auto const& matPtr : other.matrices_)
      matrices_.emplace_back(matPtr ? matPtr->clone() : nullptr);

    transition_ = other.transition_ ? other.transition_->clone() : nullptr;
  }

  // Move-semantics
  AbstractOperator(AbstractOperator&&) noexcept = default;
  AbstractOperator& operator=(const AbstractOperator&) = default;

  // Named swap for strong exception safety
  void swap(AbstractOperator& other) noexcept {
    using std::swap;
    swap(matrices_,    other.matrices_);
    swap(transition_,  other.transition_);
    swap(prevParams_,  other.prevParams_);
    swap(popIndices_,  other.popIndices_);
    // base class swap is not needed
  }
  friend void swap(AbstractOperator& a, AbstractOperator& b) noexcept {
    a.swap(b);
  }

  // Polymorphic clone
  AbstractOperator* clone() const override = 0;

  // Convenience wrapper
  std::unique_ptr<AbstractOperator> cloneOperator() const {
    return std::unique_ptr<AbstractOperator>(clone());
  }

  virtual ~AbstractOperator() {
    std::vector<std::string> paramNames(0);
    paramNames.reserve(getParameters().size());

    for(size_t i = 0; i < getParameters().size(); ++i)
      paramNames.emplace_back(getParameters()[i].getName());

    deleteParameters_(paramNames);
  }

  // Replace entire parameter set; triggers update if any changed
  void setParameters(const bpp::ParameterList& params) {
    bpp::AbstractParameterAliasable::setParametersValues(params);
  }
  void fireParameterChanged(const bpp::ParameterList& params) {
    if (matchParametersValues(params))
      updateMatrices_();
  }

  // Accessors
  const std::vector<size_t>& getPopIndices() const noexcept {
    return popIndices_;
  }
  const std::vector<std::unique_ptr<MatrixEngine>>& getMatrices() const noexcept {
    return matrices_;
  }
  const MatrixEngine& getMatrix(size_t idx) const {
    return *matrices_.at(idx);
  }

  // Extract the raw MatrixVariant — empty if no transition set
  MatrixEngine::MatrixVariant getTransitionMatrixVariant() const {
    if (!transition_) return MatrixEngine::MatrixVariant{};
    return transition_->matrixVariant();
  }

  // Extract the Eigen‐view variant (dense/sparse) — empty if no transition
  MatrixEngine::MatrixEigenVariant getTransitionMatrixVariantEigen() const {
    if (!transition_) return MatrixEngine::MatrixEigenVariant{};
    return transition_->toEigenMatrixVariant();
  }

  // Debug‐print the “delta log‐density” matrix to file
  virtual void printDeltaLDMat(const std::string& fileName);

  // Uniformly scale the assembled transition
  void scaleMatrix(double factor) {
    if (!transition_)
      throw bpp::Exception("scaleMatrix: transition matrix is not set");
    transition_->scaleMatrix(factor);
  }

protected:
  // 1) Build raw delta‐matrices (one per pop or pop‐pair)
  virtual void setUpMatrices_(const SumStatsLibrary& sslib) = 0;

  // 2) Re-scale them by current parameters
  virtual void updateMatrices_() = 0;

  // 3) Sum them into `transition_`; default implementation provided here
  virtual void assembleTransitionMatrix_() {
    if (matrices_.empty()) {
      transition_.reset();
      return;
    }

    // Clone first delta as the base
    transition_ = matrices_[0]
      ? matrices_[0]->clone()
      : nullptr;

    // Add all others via your variant‐safe helper
    for (size_t i = 1; i < matrices_.size(); ++i) {
      if (!matrices_[i]) continue;
      transition_->addToMatrix(matrices_[i]->matrixVariant());
    }
  }
};

#endif
// ==== END AbstractOperator.hpp ====
