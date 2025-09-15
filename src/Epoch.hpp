/*
 * Authors: Gustavo V. Barroso
 * Created: 30/08/2022
 * Last modified: 15/09/2025
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
#include <eigen3/Eigen/Eigenvalues>
#include <eigen3/unsupported/Eigen/MatrixFunctions>
#include <eigen3/unsupported/Eigen/MPRealSupport>

#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/AbstractParameterAliasable.h>

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

class Epoch: public bpp::AbstractParameterAliasable
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

  // for continuous-time integration
  double dt_;         // default time step for fixed integration
  size_t steps_;      // default number of steps
  double totalTime_;  // default total time for adaptive integration
  double tolerance_;  // default error tolerance

public:
  Epoch():
  bpp::AbstractParameterAliasable(""),
  name_(),
  ssl_(),
  startGen_(0),
  endGen_(0),
  pops_(0),
  operators_(0),
  engine_(std::make_unique<MatrixEngine>(false)),
  dt_(0.),
  steps_(0),
  totalTime_(0.),
  tolerance_(0.)
  { }

  Epoch(const std::string& name, const SumStatsLibrary& ssl, size_t start, size_t end,
        const std::vector<std::shared_ptr<Population>>& pops,
        const std::vector<std::shared_ptr<AbstractOperator>>& ops,
        double dt = 1e-3, size_t steps = 1e+3, double totalTime = 1., double tol = 1e-6):
  bpp::AbstractParameterAliasable(""),
  name_(name),
  ssl_(ssl),
  startGen_(start),
  endGen_(end),
  pops_(pops),
  operators_(ops),
  engine_(std::make_unique<MatrixEngine>(false)),
  dt_(dt),
  steps_(steps),
  totalTime_(totalTime),
  tolerance_(tol)
  {
    for(auto it = std::begin(operators_); it != std::end(operators_); ++it)
      addParameters_((*it)->getParameters());

    bpp::AbstractParameterAliasable::setNamespace(name + ".");
    init_();
  }

  ~Epoch()
  {
    std::vector<std::string> paramNames(0);
    paramNames.reserve(getParameters().size());

    for(size_t i = 0; i < getParameters().size(); ++i)
      paramNames.emplace_back(getParameters()[i].getName());

    deleteParameters_(paramNames);
  }

  Epoch* clone() const
  {
    return new Epoch(*this);
  }

  void fireParameterChanged(const bpp::ParameterList& params);

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

  const std::string& getName()
  {
    return name_;
  }

  size_t start()
  {
    return startGen_;
  }

  size_t end()
  {
    return endGen_;
  }

  size_t duration()
  {
    return startGen_ - endGen_;
  }

  double getDt()
  {
    return dt_;
  }

  size_t getSteps()
  {
    return steps_;
  }

  double getTotalTime()
  {
    return totalTime_;
  }

  double getTolerance()
  {
    return tolerance_;
  }

  auto getTransitionMatrix() const -> MatrixEngine::SparseMatrixVariant
  {
    return engine_->getRawMatrix();
  }

  auto getSteadyStateVector() const -> MatrixEngine::VectorVariantEigen
  {
    return engine_->getRawVector();
  }

  size_t getNumPops()
  {
    return pops_.size();
  }

  void setDt(double dt)
  {
    dt_ = dt;
  }

  void setNumSteps(size_t numSteps)
  {
    steps_ = steps;
  }

  void setTotalTime(double time)
  {
    totalTime_ = time;
  }

  void setTolerance(double tol)
  {
    tolerance_ = tol;
  }

  const std::vector<std::shared_ptr<Population>>& getPops()
  {
    return pops_;
  }

  void printAttributes(std::ostream& stream)
  {
    stream << name_ << ", from " << startGen_ << " to " << endGen_ << "\n";

    for(auto it = std::begin(pops_); it != std::end(pops_); ++it)
    {
      stream << "\t";
      (*it)->printAttributes(stream);
    }
  }

  std::shared_ptr<Population> fetchPop(size_t id)
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

  std::shared_ptr<Population> fetchPop(const std::string& name)
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

  void computeExpectedSumStatsDiscrete(const MatrixEngine::VectorVariantEigen& y);

  void transferStatistics(MatrixEngine::VectorVariantEigen& y) const;

  void updateMoments(const MatrixEngine::VectorVariantEigen& y);

  void printMoments(std::ostream& stream);

  void printMomentsIntermediate(
  MatrixEngine::VectorVariantEigen& y,
  const std::string& modelName,
  size_t interval,
  const std::vector<std::string>& momNames);

  void printRecursions(std::ostream& stream);

  void printTransitionMat(const std::string& fileName) const;

  void computePseudoSteadyStateDiscrete();

  void computeEigenSteadyState();

  void testSteadyState();

  template<typename Scalar>
  Eigen::Matrix<Scalar, Eigen::Dynamic, 1> integrateTyped(const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& moms) const;

  MatrixEngine::VectorVariantEigen integrate(const MatrixEngine::VectorVariantEigen& moms) const;

  template<typename Scalar>
  Eigen::Matrix<Scalar, Eigen::Dynamic, 1> integrateAdaptiveTyped(
  const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& moms,
  double dtMin = 1e-6,
  double dtMax = 1.0) const;

  MatrixEngine::VectorVariantEigen integrateAdaptive(const MatrixEngine::VectorVariantEigen& moms) const;

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
  double fetchNu(size_t popId)
  {
    return fetchNu(popId, pops_[popId]->getParent()->getSize());
  }

  inline double fetchConditionNumber() const
  {
    return std::visit([](const auto& mat)
    {
      using Scalar = typename std::decay_t<decltype(mat)>::Scalar;
      using DenseMatrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

      DenseMatrix dense(mat.mat_); // convert sparse to dense
      Eigen::JacobiSVD<DenseMatrix> svd(dense);

      const auto& singularValues = svd.singularValues();
      return singularValues(0).toDouble() / singularValues(singularValues.size() - 1).toDouble();
    }, engine_->getMatrixVariant());
  }

  inline EigenResult findLeadingEigenpair() const
  {
    return std::visit([](const auto& mat) -> EigenResult
    {
      using Scalar = typename std::decay_t<decltype(mat)>::Scalar;
      using DenseMatrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

      DenseMatrix dense(mat.mat_); // convert sparse to dense
      Eigen::EigenSolver<DenseMatrix> es(dense);

      // Find index of largest real eigenvalue
      size_t idx = 0;
      for(size_t i = 1; i < es.eigenvalues().size(); ++i)
      {
        if(es.eigenvalues().real()(i) > es.eigenvalues().real()(idx))
          idx = i;
      }

      // Extract and normalize the corresponding eigenvector
      auto vec = es.eigenvectors().col(idx).real();
      vec.normalize();

      // Convert to mpreal for consistency
      Eigen::Matrix<mpfr::mpreal, Eigen::Dynamic, 1> vecMP(vec.size());
      for (size_t i = 0; i < vec.size(); ++i)
        vecMP(i) = mpfr::mpreal(vec(i));

      return { idx, mpfr::mpreal(es.eigenvalues().real()(idx)), vecMP };
    }, engine_->getMatrixVariant());
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
