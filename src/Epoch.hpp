/*
 * Authors: Gustavo V. Barroso
 * Created: 30/08/2022
 * Last modified: 21/10/2025
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

#define EIGEN_DONT_PARALLELIZE 0
#define EIGEN_USE_THREADS

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/SVD>
#include <eigen3/Eigen/Eigenvalues>
#include <eigen3/unsupported/Eigen/MatrixFunctions>
#include <eigen3/unsupported/Eigen/MPRealSupport>
#include <eigen3/unsupported/Eigen/CXX11/ThreadPool>

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wnon-virtual-dtor"
#include <Spectra/SymEigsSolver.h>
#include <Spectra/GenEigsSolver.h>
#include <Spectra/MatOp/DenseGenMatProd.h>
#pragma clang diagnostic pop

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
        bool multiplyOperators, // whether to obtain transition matrix by multiplying (or summing) operator matrices
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
    init(multiplyOperators);
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
    // init may be optional if engine_ already cloned appropriately
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
                                const std::string& modelName,
                                size_t interval,
                                const std::vector<std::string>& momNames);

  void printRecursions(std::ostream& stream);

  void printTransitionMat(const std::string& fileName) const;

  void computePowerSteadyStateDiscrete(double tol = 1e-16);

  void computePowerSteadyStateContinuous(double burnInTime = 0.1, double dt = 1e-3, double tol = 1e-16);

  // typed helper for the continuous power steady-state solver
  template<typename Scalar>
  bool computePowerSSContinuousTyped(
    const Matrix<Scalar>& A,
    double burnInTime,
    double dt,
    double tol,
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& outY) const;

  template<typename Scalar>
  bool computePowerSSAdaptiveTyped(
    const Matrix<Scalar>& A,
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
    std::cout << "\nCondition Number for Transition Matrix, epoch " << name_ << " = " << cond << "\n";
  }

  // computes and returns relative population size (use for continuous-time integration)
  double fetchNu(size_t popId, double Nref)
  {
    double Nfocal = pops_[popId]->getSize();
    return Nfocal / Nref;
  }

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

  // finds leading Eigen pair using Eigen 3.4 routines
  inline EigenResult findLeadingEigenpair() const
  {
    auto const& eigenVar = engine_->toEigenMatrixVariant();

    return std::visit(overloaded {[this](auto const& M) -> EigenResult
    {
      using DenseD = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
      using VecD = Eigen::Matrix<double, Eigen::Dynamic, 1>;

      // converts to dense<double>
      DenseD dense = M.template cast<double>();

      // solve eigenproblem
      Eigen::EigenSolver<DenseD> es(dense);
      if(es.info() != Eigen::Success)
        throw bpp::Exception("findLeadingEigenpair(): EigenSolver failed");

      // picks largest real eigenvalue
      auto evals = es.eigenvalues().real();
      Eigen::Index idx = 0;

      for(Eigen::Index i = 1; i < evals.size(); ++i)
      {
        if(evals(i) > evals(idx))
          idx = i;
      }

      std::cout << std::setprecision(12) << "\nleading eigenval (base Eigen 3.4) =  " << evals(idx) << "\n";

      if(es.eigenvalues().real()(idx) > 1. + 1e-5 || es.eigenvalues().real()(idx) < 1. - 1e-5)
      {
        double cond = fetchConditionNumber();
        std::cout << "\nCondition Number of transition matrix = " << cond << "\n";
        throw bpp::Exception("Epoch::Bad Leading Eigenvalue! Consider using a smaller order of 1-2p factors.\n");
      }

      VecD vecD = es.eigenvectors().col(idx).real();
      // I moment embodies scaling constant used by Eigen
      vecD /= vecD(ssl_.findCompressedIndex(ssl_.getMoment("I")));

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

  // finds leading Eigen pair using Spectra as a multi-threaded backend
  // this improves computationally considerably because we only want the leading eigenpar
  // (don't need full Eigen3.4 decomposition)
  inline EigenResult findLeadingEigenpairSpectra() const
  {
    auto const& eigenVar = engine_->toEigenMatrixVariant();

    return std::visit(overloaded {[this](auto const& M) -> EigenResult
    {
      using DenseD = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
      using VecD = Eigen::Matrix<double, Eigen::Dynamic, 1>;

      // converts to dense<double>
      DenseD dense = M.template cast<double>();

      // wraps matrix in Spectra operator
      Spectra::DenseGenMatProd<double> op(dense);

      // creates solver: compute 1 largest eigenvalue
      Spectra::GenEigsSolver<Spectra::DenseGenMatProd<double>> eigs(op, 1, std::min<int>(20, dense.cols()));

      size_t pop = ssl_.getPopIndices()[0];
      double mu = getParameterValue("u_" + bpp::TextTools::toString(pop));
      double s = getParameterValue("s_" + bpp::TextTools::toString(pop));
      size_t twoN = static_cast<size_t>(1. / getParameterValue("1/2N_" + bpp::TextTools::toString(pop)));

      double hr = twoN * mu;
      double hl = (std::abs(s) < 1e-14) ? twoN * mu : (2. * twoN * mu * std::exp(2. * twoN * s) / (std::exp(2. * twoN * s) - 1.0) - mu / s);

      VecD y(dense.cols()); // initial guess

      double f = 1.0; // scaling factor to mimic the factors of 1-2p
      for(Eigen::Index i = 0; i < y.size(); ++i)
      {
        const auto& prefix = ssl_.getBasis()[size_t(i)]->getPrefix();

        if(i > 0 && prefix == ssl_.getBasis()[size_t(i - 1)]->getPrefix())
          f *= 0.925; // same kind of moment, applies heuristic decay

        else
          f = 1.0; // resets for next moment

        if(prefix == "Hl")
          y(i) = hl * f;

        else if(prefix == "Hr")
          y(i) = hr * f;

        else if(prefix == "pi2")
          y(i) = hr * hl * f * 1.3;

        else if(prefix == "I")
          y(i) = 1.0;

        else if(prefix == "DD")
          y(i) = hr * hl * f * 5e-1;

        else // if(prefix == "Dr")
          y(i) = hr * hl * f * 3e-1;
      }

      eigs.init(y.data()); // helps convergence
      int nconv = eigs.compute(); // default = largest magnitude eigenvlaue

      if(nconv < 1 || eigs.info() != Spectra::CompInfo::Successful)
        throw bpp::Exception("findLeadingEigenpair(): Spectra solver failed");

      // extracts eigenvalue and eigenvector
      Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> vecC = eigs.eigenvectors().col(0);
      VecD vecD(vecC.size());

      for(Eigen::Index i = 0; i < vecC.size(); ++i)
        vecD(i) = vecC(i).real();

      // I moment embodies scaling constant used by Eigen
      vecD /= vecD(ssl_.findCompressedIndex(ssl_.getMoment("I")));

      double lambda = eigs.eigenvalues()(0).real(); // handles complex return

      std::cout << std::setprecision(12) << "\nleading eigenval (Spectra) =  " << lambda << "\n";

      if(lambda > 1. + 1e-5 || lambda < 1. - 1e-5)
      {
        double cond = fetchConditionNumber();
        std::cout << "\nCondition Number of transition matrix = " << cond << "\n";
        throw bpp::Exception("Epoch::Bad Leading Eigenvalue! Consider using a smaller order of 1-2p factors.\n");
      }

      // converts to high precision
      Eigen::Matrix<mpfr::mpreal, Eigen::Dynamic, 1> vecMP(vecD.size());
      for(Eigen::Index i = 0; i < vecD.size(); ++i)
        vecMP(i) = mpfr::mpreal(vecD(i));

      return EigenResult
      {
        0, // Spectra returns one eigenpair, index is always 0
        mpfr::mpreal(lambda),
        vecMP
      };

    }
    }, eigenVar);
  }

  void init(bool multiplyOperators);

private:
  void updateOperators_(const bpp::ParameterList& params)
  {
    for(auto it = std::begin(operators_); it != std::end(operators_); ++it)
      (*it)->fireParameterChanged(params);
  }
};

#endif
