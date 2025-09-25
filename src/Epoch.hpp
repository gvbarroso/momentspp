/*
 * Authors: Gustavo V. Barroso
 * Created: 30/08/2022
 * Last modified: 23/09/2025
 *
 */

#ifndef _EPOCH_H_
#define _EPOCH_H_

#include <iostream>
#include <vector>
#include <memory>
#include <utility>
#include <algorithm>
#include <map>

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/SVD>
#include <eigen3/Eigen/Eigenvalues>
#include <eigen3/unsupported/Eigen/MatrixFunctions>
#include <eigen3/unsupported/Eigen/MPRealSupport>

#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/AbstractParameterAliasable.h>

#include "VariantUtils.hpp"
#include "EpochIntegrators.hpp"
#include "AbstractOperator.hpp"
#include "Admixture.hpp"
#include "Mutation.hpp"
#include "SumStatsLibrary.hpp"
#include "Population.hpp"
#include "MatrixEngine.hpp"

struct EigenResult
{
  size_t index;
  mpfr::mpreal value;
  Eigen::Matrix<mpfr::mpreal, Eigen::Dynamic, 1> vector;
};

class Epoch : public bpp::AbstractParameterAliasable
{

private:
  std::string name_;
  SumStatsLibrary ssl_; // *this epoch has its own set of moments using its population indices

  // generations ago, from past to present (startGen_ > endGen_)
  size_t startGen_;
  size_t endGen_;

  std::vector<std::shared_ptr<Population>> pops_;
  std::vector<std::shared_ptr<AbstractOperator>> operators_; // each operator contains matrices and a subset of the parameters

  // engine_ holds the steady state vector as well as all sparse operators summed into a Sparse matrix
  std::unique_ptr<MatrixEngine> engine_;

public:
  Epoch():
  bpp::AbstractParameterAliasable(""),
  name_(),
  ssl_(),
  startGen_(0),
  endGen_(0),
  pops_(0),
  operators_(0),
  engine_(std::make_unique<MatrixEngine>(false))
  { }

  Epoch(const std::string& name, const SumStatsLibrary& ssl, size_t start, size_t end,
        const std::vector<std::shared_ptr<Population>>& pops,
        const std::vector<std::shared_ptr<AbstractOperator>>& ops):
  bpp::AbstractParameterAliasable(""),
  name_(name),
  ssl_(ssl),
  startGen_(start),
  endGen_(end),
  pops_(pops),
  operators_(ops),
  engine_(std::make_unique<MatrixEngine>(false))
  {
    for(auto it = std::begin(operators_); it != std::end(operators_); ++it)
      addParameters_((*it)->getParameters());

    bpp::AbstractParameterAliasable::setNamespace(name + ".");
    init_();
  }

  Epoch(const Epoch& other):
  bpp::AbstractParameterAliasable(""),
  name_(other.name_),
  ssl_(other.ssl_),
  startGen_(other.startGen_),
  endGen_(other.endGen_),
  pops_(other.pops_),
  operators_(),
  engine_(other.engine_ ? other.engine_->clone() : nullptr)
  {
    operators_.reserve(other.operators_.size());
    for(const auto& op : other.operators_)
      operators_.push_back(op ? std::shared_ptr<AbstractOperator>(op->clone()) : nullptr);

    // copy parameters and namespace exactly as constructor does
    for(const auto& op : operators_)
    {
      if(op)
        addParameters_(op->getParameters());
    }

    bpp::AbstractParameterAliasable::setNamespace(name_ + ".");
    // init_ may be optional if engine_ already cloned appropriately
  }

  Epoch& operator=(const Epoch& other)
  {
    if(this != &other)
    {
      // clear current parameters registered under this namespace
      std::vector<std::string> paramNames;
      paramNames.reserve(getParameters().size());

      for(size_t i = 0; i < getParameters().size(); ++i)
        paramNames.emplace_back(getParameters()[i].getName());

      deleteParameters_(paramNames);

      name_ = other.name_;
      ssl_ = other.ssl_;
      startGen_ = other.startGen_;
      endGen_ = other.endGen_;
      pops_ = other.pops_;

      // clone operators
      operators_.clear();
      operators_.reserve(other.operators_.size());
      for(const auto& op : other.operators_)
      {
        if(op)
          operators_.push_back(std::shared_ptr<AbstractOperator>(op->clone()));
        else
          operators_.push_back(nullptr);
      }

      // clone engine
      engine_ = other.engine_ ? other.engine_->clone() : nullptr;

      // re-register parameters
      for(const auto& op : operators_)
      {
        if(op)
          addParameters_(op->getParameters());
      }

      bpp::AbstractParameterAliasable::setNamespace(name_ + ".");
    }

    return *this;
  }

  Epoch(Epoch&&) noexcept = default;
  Epoch& operator=(Epoch&&) noexcept = default;

  ~Epoch()
  {
    std::vector<std::string> paramNames(0);
    paramNames.reserve(getParameters().size());

    for(size_t i = 0; i < getParameters().size(); ++i)
      paramNames.emplace_back(getParameters()[i].getName());

    deleteParameters_(paramNames);
  }

  Epoch* clone() const override
  {
    return new Epoch(*this);
  }

  void fireParameterChanged(const bpp::ParameterList& params) override;

  void setParameters(const bpp::ParameterList& params)
  {
    bpp::AbstractParameterAliasable::setParametersValues(params);
  }

  const MatrixEngine& getEngine() const
  {
    return *engine_;
  }

  MatrixEngine& getEngine()
  {
    return *engine_;
  }

  const std::string& getName() const
  {
    return name_;
  }

  size_t start() const
  {
    return startGen_;
  }

  size_t end() const
  {
    return endGen_;
  }

  size_t duration() const
  {
    return startGen_ - endGen_;
  }

  auto getTransitionMatrix() const -> MatrixEngine::MatrixEigenVariant
  {
    return engine_->toEigenMatrixVariant();
  }

  auto getSteadyStateVector() const -> MatrixEngine::VectorVariant
  {
    return engine_->getVectorVariant();
  }

  size_t getNumPops() const
  {
    return pops_.size();
  }

  const std::vector<std::shared_ptr<Population>>& getPops()
  {
    return pops_;
  }

  void printAttributes(std::ostream& stream)
  {
    stream << name_ << ", from " << startGen_ << " to " << endGen_ << "\n";

    for(const auto& pop : pops_)
    {
      stream << "\t";
      pop->printAttributes(stream);
    }
  }

  std::shared_ptr<Population> fetchPop(size_t id) const
  {
    std::shared_ptr<Population> pop = nullptr;
    for(auto it = std::begin(pops_); it != std::end(pops_); ++it)
    {
      if((*it)->getId() == id)
        pop = (*it);
    }

    assert(pop != nullptr);
    return pop;
  }

  std::shared_ptr<Population> fetchPop(const std::string& name) const
  {
    std::shared_ptr<Population> pop = nullptr;
    for(auto it = std::begin(pops_); it != std::end(pops_); ++it)
    {
      if((*it)->getName() == name)
        pop = (*it);
    }

    assert(pop != nullptr);
    return pop;
  }

  const SumStatsLibrary& getSslib() const
  {
    return ssl_;
  }

  SumStatsLibrary& getSslib()
  {
    return ssl_;
  }

  const std::vector<std::shared_ptr<Moment>>& getMoments() const
  {
    return ssl_.getMoments();
  }

  std::vector<std::shared_ptr<Moment>>& getMoments()
  {
    return ssl_.getMoments();
  }

  const std::vector<std::shared_ptr<Moment>>& getBasis() const
  {
    return ssl_.getBasis();
  }

  std::vector<std::shared_ptr<Moment>>& getBasis()
  {
    return ssl_.getBasis();
  }

  std::vector<size_t> fetchSelectedPopIds(); // for *this epoch

  void computeExpectedSumStatsDiscrete(MatrixEngine::VectorVariant& y);

  void transferStatistics(MatrixEngine::VectorVariant& y);

  void updateMoments(const MatrixEngine::VectorVariant& y);

  void printMoments(std::ostream& stream);

  void printMomentsIntermediate(MatrixEngine::VectorVariant& y,
                                const std::string& modelName, size_t interval,
                                const std::vector<std::string>& momNames);

  void printRecursions(std::ostream& stream);

  void printTransitionMat(const std::string& fileName) const;

  void computePseudoSteadyStateDiscrete(double tol = 1e-6);

  void computePseudoSteadyStateContinuous(double burnInTime = 0.1, double dt = 1e-3, double tol = 1e-6);

  // typed helper for the continuous pseudo‐steady solver
  template<typename Scalar>
  bool computePseudoSSContinuousTyped(
    const Eigen::SparseMatrix<Scalar>& A,
    double burnInTime,
    double dt,
    double tol,
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& outY) const;

  void computeEigenSteadyState();

  void testSteadyState();

  MatrixEngine::VectorVariant integrate(double dt, double totalTime) const;

  MatrixEngine::VectorVariant integrateAdaptive(double dt, double totalTime, double tol, double dtMin = 1e-6, double dtMax = 1e-2) const;

  void printConditionNumber()
  {
    double cond = fetchConditionNumber();
    std::cout << "Condition Number for Transition Matrix, epoch " << name_ << " = " << cond << "\n";
  }

  // computes and returns relative population size (use for continuous-time integration)
  double fetchNu(size_t popId, double Nref)
  {
    double Nfocal = pops_[popId]->getSize();
    return Nfocal / Nref;
  }

  // computes and returns relative population size w.r.t "same pop" in previous epoch
  /*double fetchNu(size_t popId)
  {
    return fetchNu(popId, pops_[popId]->getParent()->getSize());
  }*/

  inline double fetchConditionNumber() const
  {
    auto const& eigenVar = engine_->toEigenMatrixVariant();

    return std::visit(overloaded {[](auto const& M) -> double
    {
      using Dense = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;

      // sparse<double> or sparse<mpreal> → sparse<double> → dense<double>
      Dense dense = M.template cast<double>();
      Eigen::JacobiSVD<Dense> svd(dense, Eigen::ComputeThinU | Eigen::ComputeThinV);

      if(svd.info() != Eigen::Success)
        throw bpp::Exception("fetchConditionNumber(): SVD failed");

      auto s = svd.singularValues();
      return (s.size() > 1) ? static_cast<double>(s(0) / s(s.size() - 1)) : 0.0;
    }
  }, eigenVar);
}

  inline EigenResult findLeadingEigenpair() const
  {
    auto const& eigenVar = engine_->toEigenMatrixVariant();

    return std::visit(overloaded {[](auto const& M) -> EigenResult
    {
      using DenseD = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
      using VecD   = Eigen::Matrix<double, Eigen::Dynamic, 1>;

      // to dense<double>
      DenseD dense = M.template cast<double>();

      // solve eigenproblem
      Eigen::EigenSolver<DenseD> es(dense);
      if (es.info() != Eigen::Success)
        throw bpp::Exception("findLeadingEigenpair(): EigenSolver failed");

      // pick largest real eigenvalue
      auto evals = es.eigenvalues().real();
      Eigen::Index idx = 0;

      for(Eigen::Index i = 1; i < evals.size(); ++i)
      {
        if(evals(i) > evals(idx))
          idx = i;
      }

      // normalize its eigenvector
      VecD vecD = es.eigenvectors().col(idx).real();
      vecD.normalize();

      // convert to high precision
      Eigen::Matrix<mpfr::mpreal, Eigen::Dynamic, 1> vecMP(vecD.size());
      for(Eigen::Index i = 0; i < vecD.size(); ++i)
        vecMP(i) = mpfr::mpreal(vecD(i));

      return EigenResult
      {
        size_t(idx),
        mpfr::mpreal(evals(idx)),
        vecMP
      };
    }
  }, eigenVar);
}

private:
  void init_();

  void updateOperators_(const bpp::ParameterList& params)
  {
    for(auto it = std::begin(operators_); it != std::end(operators_); ++it)
      (*it)->fireParameterChanged(params);
  }
};

#endif
