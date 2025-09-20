#include "Epoch.hpp"
#include "EpochIntegrators.hpp"

using bpp::Exception;

//------------------------------------------------------------------------------
// Construction
//------------------------------------------------------------------------------

Epoch::Epoch(std::string name,
             size_t startGen,
             std::unique_ptr<MatrixEngine> engine,
             std::vector<std::shared_ptr<AbstractOperator>> operators,
             SumStatsLibrary ssl)
  : name_(std::move(name)),
    startGen_(startGen),
    engine_(std::move(engine)),
    operators_(std::move(operators)),
    ssl_(std::move(ssl))
{
  init_();
}

//------------------------------------------------------------------------------
// init_: accumulate transition matrices & prep engine
//------------------------------------------------------------------------------

void Epoch::init_()
{
  if (operators_.empty())
    throw Exception("Epoch::init_() called with no operators.");

  if (!engine_)
    throw Exception("Epoch::init_() called with null engine_.");

  // start with first operator’s matrix
  auto acc = operators_.front()->getTransitionMatrixVariantEigen();

  // accumulate only matching scalar types
  for (size_t i = 1; i < operators_.size(); ++i) {
    auto next = operators_[i]->getTransitionMatrixVariantEigen();
    visit_same_scalar(acc, next, [&](auto& A, auto const& B) {
      A += B;
    });
  }

  engine_->setMatrixFromEigen(std::move(acc));
  engine_->addIdentityInPlace();
  engine_->pruneInPlace();
  engine_->compressInPlace();
}

//------------------------------------------------------------------------------
// fireParameterChanged: update operators & rebuild matrix
//------------------------------------------------------------------------------

void Epoch::fireParameterChanged(const bpp::ParameterList& params)
{
  for (auto& op : operators_)
    op->fireParameterChanged(params);

  init_();
}

//------------------------------------------------------------------------------
// Steady‐state via leading eigenpair
//------------------------------------------------------------------------------

double Epoch::fetchConditionNumber() const
{
  auto matVar = engine_->toEigenMatrixVariant();
  return std::visit([](auto const& M) -> double {
    using Dense = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
    Dense dense = M.template cast<double>();
    Eigen::JacobiSVD<Dense> svd(dense,
      Eigen::ComputeThinU | Eigen::ComputeThinV);
    if (svd.info() != Eigen::Success)
      throw Exception("fetchConditionNumber(): SVD failed");
    auto s = svd.singularValues();
    return (s.size()>1) ? s(0)/s(s.size()-1) : 0.0;
  }, matVar);
}

EigenResult Epoch::findLeadingEigenpair() const
{
  auto matVar = engine_->toEigenMatrixVariant();
  return std::visit([](auto const& M) -> EigenResult {
    using DenseD = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
    using VecD   = Eigen::Matrix<double, Eigen::Dynamic, 1>;
    DenseD dense = M.template cast<double>();
    Eigen::EigenSolver<DenseD> es(dense);
    if (es.info() != Eigen::Success)
      throw Exception("findLeadingEigenpair(): EigenSolver failed");
    auto evals = es.eigenvalues().real();
    Eigen::Index idx = 0;
    for (Eigen::Index i = 1; i < evals.size(); ++i)
      if (evals(i) > evals(idx)) idx = i;
    VecD vecD = es.eigenvectors().col(idx).real();
    vecD.normalize();
    Eigen::Matrix<mpfr::mpreal, Eigen::Dynamic, 1> vecMP(vecD.size());
    for (Eigen::Index i = 0; i < vecD.size(); ++i)
      vecMP(i) = mpfr::mpreal(vecD(i));
    return EigenResult{
      size_t(idx),
      mpfr::mpreal(evals(idx)),
      vecMP
    };
  }, matVar);
}

void Epoch::computeEigenSteadyState()
{
  EigenResult res = findLeadingEigenpair();
  if (res.value.toDouble() > 1.0 + 1e-14)
    throw Exception("Leading eigenvalue > 1");

  auto vecVar = engine_->toEigenVectorVariant();
  std::visit([&](auto& dst) {
    using Scalar = typename std::decay_t<decltype(dst)>::Scalar;
    for (Eigen::Index i = 0; i < dst.size(); ++i)
      dst(i) = static_cast<Scalar>(res.vector(i));
  }, vecVar);

  engine_->setVectorFromEigen(std::move(vecVar));
  updateMoments(engine_->toEigenVectorVariant());
}

//------------------------------------------------------------------------------
// Discrete‐time pseudo‐steady‐state (power method)
//------------------------------------------------------------------------------

void Epoch::computePseudoSteadyStateDiscrete(double tol)
{
  auto matVar = engine_->toEigenMatrixVariant();
  bool converged = false;

  std::visit([&](auto const& A) {
    using Scalar = typename std::decay_t<decltype(A)>::Scalar;
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> y(A.cols());
    y.setOnes();
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> prev;

    for (size_t iter = 0; iter < 10000; ++iter) {
      prev = y;
      y = A * y;
      double maxRel = 0;
      for (Eigen::Index i = 0; i < y.size(); ++i)
        maxRel = std::max(maxRel,
          std::abs(static_cast<double>(y(i)-prev(i)))
          / std::max(1.0, std::abs(static_cast<double>(prev(i)))));
      if (maxRel <= tol) { converged = true; break; }
    }

    if (converged) {
      auto vecVar = engine_->toEigenVectorVariant();
      std::visit([&](auto& dst) {
        using DScalar = typename std::decay_t<decltype(dst)>::Scalar;
        if constexpr (std::is_same_v<DScalar, Scalar>)
          dst = y;
      }, vecVar);
      engine_->setVectorFromEigen(std::move(vecVar));
    }
  }, matVar);

  if (!converged)
    computeEigenSteadyState();

  updateMoments(engine_->toEigenVectorVariant());
}

//------------------------------------------------------------------------------
// Continuous‐time pseudo‐steady‐state (Crank–Nicolson)
//------------------------------------------------------------------------------

template<typename Scalar>
bool Epoch::computePseudoSSContinuousTyped(
    const Eigen::SparseMatrix<Scalar>& A,
    double burnInTime,
    double dt,
    double tol,
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& outY) const
{
  Eigen::Matrix<Scalar, Eigen::Dynamic, 1> y(A.cols());
  y.setOnes();

  size_t burnSteps = size_t(burnInTime / dt);
  for (size_t i = 0; i < burnSteps; ++i)
    y = A * y;

  Eigen::Matrix<Scalar, Eigen::Dynamic, 1> prev;
  for (size_t iter = 0; iter < burnSteps*10; ++iter) {
    prev = y;
    y = A * y;
    double maxRel = 0;
    for (Eigen::Index i = 0; i < y.size(); ++i)
      maxRel = std::max(maxRel,
        std::abs(static_cast<double>(y(i)-prev(i)))
        / std::max(1.0, std::abs(static_cast<double>(prev(i)))));
    if (maxRel <= tol) { outY = y; return true; }
  }
  return false;
}

void Epoch::computePseudoSteadyStateContinuous(double burnInTime,
                                               double dt,
                                               double tol)
{
  auto matVar = engine_->toEigenMatrixVariant();
  bool converged = false;
  Eigen::Matrix<double, Eigen::Dynamic, 1> yD;
  Eigen::Matrix<mpfr::mpreal, Eigen::Dynamic, 1> yM;

  std::visit([&](auto const& A) {
    using Scalar = typename std::decay_t<decltype(A)>::Scalar;
    if constexpr (std::is_same_v<Scalar,double>)
      converged = computePseudoSSContinuousTyped(A, burnInTime, dt, tol, yD);
    else
      converged = computePseudoSSContinuousTyped(A, burnInTime, dt, tol, yM);
  }, matVar);

  if (!converged) {
    computeEigenSteadyState();
    return;
  }

  MatrixEngine::VectorVariantEigen vecVar =
    (yD.size() ? MatrixEngine::VectorVariantEigen(std::move(yD))
               : MatrixEngine::VectorVariantEigen(std::move(yM)));

  engine_->setVectorFromEigen(std::move(vecVar));
  updateMoments(engine_->toEigenVectorVariant());
}

//------------------------------------------------------------------------------
// Integration wrappers
//------------------------------------------------------------------------------

MatrixEngine::VectorVariant Epoch::integrate(double dt, double totalTime) const
{
  const auto& matVar = engine_->getMatrixVariant();
  const auto& vecVar = engine_->getVectorVariant();
  MatrixEngine::VectorVariant result;

  visit_same_scalar(matVar, vecVar, [&](auto const& M, auto const& V) {
    using T = typename std::decay_t<decltype(M)>::Scalar;
    if constexpr (std::is_same_v<T,double>)
      result = integrateDoubleCN(M, V, dt, totalTime);
    else
      result = integrateMpfrCN(M, V, dt, totalTime);
  });

  return result;
}

MatrixEngine::VectorVariant Epoch::integrateAdaptive(double dt,
                                                     double totalTime,
                                                     double tol,
                                                     double dtMin,
                                                     double dtMax) const
{
  const auto& matVar = engine_->getMatrixVariant();
  const auto& vecVar = engine_->getVectorVariant();
  MatrixEngine::VectorVariant result;

  visit_same_scalar(matVar, vecVar, [&](auto const& M, auto const& V) {
    using T = typename std::decay_t<decltype(M)>::Scalar;
    if constexpr (std::is_same_v<T,double>)
      result = integrateAdaptiveDoubleCN(M, V, dt, totalTime, tol, dtMin, dtMax);
    else
      result = integrateAdaptiveMpfrCN(M, V, dt, totalTime, tol, dtMin, dtMax);
  });

  return result;
}

//------------------------------------------------------------------------------
// Moments & diagnostics
//------------------------------------------------------------------------------

void Epoch::updateMoments(const VectorVariantEigen& y)
{
  std::visit([&](auto const& vec) {
    for (size_t i = 0; i < ssl_.getBasis().size(); ++i)
      ssl_.getBasis()[i]->setValue(
        static_cast<double>(vec(Eigen::Index(i))));
  }, y);
}

void Epoch::transferStatistics(VectorVariantEigen& y) const
{
  auto vecVar = engine_->toEigenVectorVariant();
  visit_same_scalar(vecVar, y, [&](auto const& src, auto& dst) {
    if (src.size() != dst.size())
      throw Exception("Epoch::transferStatistics: size mismatch");
    dst = src;
  });
  engine_->setVectorFromEigen(std::move(vecVar));
}

void Epoch::printConditionNumber() const
{
  std::cout << "Condition Number for epoch " << name_
            << " = " << fetchConditionNumber() << "\n";
}

void Epoch::printRecursions(std::ostream& stream) const
{
  stream << "\n";
  const auto& basis = ssl_.getBasis();
  for (size_t i = 0; i < basis.size(); ++i) {
    if (basis[i]->getName()=="I") continue;
    Eigen::Index pos = basis[i]->getPosition();
    stream << "\u0394[" << basis[i]->getName() << "] = ";
    bool first = true;
    for (auto& op : operators_) {
      auto params = op->getParameters();
      for (size_t k = 0; k < params.size(); ++k) {
        double val = params[k].getValue();
        if (val==0.0) continue;
        auto wrapVar = op->getMatrix(k).getMatrixVariant();
        MatrixVariantEigen matVarEigen = std::visit(
          [](auto const& Mw){ return Mw.eigen(); }, wrapVar);
        std::visit([&](auto const& M){
          using Scalar = typename std::decay_t<decltype(M)>::Scalar;
          for (Eigen::Index l=0; l<M.cols(); ++l) {
            Scalar c = M.coeff(pos,l) / Scalar(val);
            if (c==Scalar(0)) continue;
            if (!first && c>Scalar(0)) stream<<"+";
            stream << std::fixed<<std::setprecision(3)
                   << c<<"*"<<params[k].getName()<<"*"
                   << basis[static_cast<size_t>(l)]->getName()<<" ";
            first = false;
          }
        }, matVarEigen);
      }
    }
    stream<<"\n";
  }
}

void Epoch::printTransitionMat(const std::string& fileName) const
{
  auto matVar = engine_->toEigenMatrixVariant();
  std::visit([&](auto const& M){
    std::ofstream out(fileName, std::ios::binary);
    if (!out) throw Exception("printTransitionMat: cannot open "+fileName);
    auto oldF = out.flags(); auto oldP = out.precision();
    out<<std::scientific<<std::setprecision(16)
       <<M.rows()<<" "<<M.cols()<<" "<<M.nonZeros()<<"\n";
    for (Eigen::Index k=0; k<M.outerSize(); ++k)
      for (typename std::decay_t<decltype(M)>::InnerIterator it(M,k); it; ++it)
        out<<it.row()<<" "<<it.col()<<" "<<it.value()<<"\n";
    out.flags(oldF); out.precision(oldP);
  }, matVar);
}

void Epoch::printMoments(std::ostream& stream) const
{
  for (auto& m : ssl_.getBasis())
    stream<<std::setprecision(16)
          <<m->getName()<<" = "<<m->getValue()<<"\n";
}

void Epoch::printMomentsIntermediate(VectorVariantEigen& y,
                                     const std::string& modelName,
                                     size_t interval,
                                     const std::vector<std::string>& momNames)
{
  transferStatistics(y);
  std::string fileName = modelName + "_" + name_ + "_moments.csv.gz";
  std::ofstream raw(fileName, std::ios::binary);
  boost::iostreams::filtering_ostream fout;
  fout.push(boost::iostreams::gzip_compressor());
  fout.push(raw);

  const auto& basis = ssl_.getBasis();
  size_t steps = startGen_/interval + 1;
  std::vector<size_t> indices;
  for (size_t i=0;i<basis.size();++i)
    if (std::find(momNames.begin(),momNames.end(),basis[i]->getName())
        != momNames.end())
      indices.push_back(i);

  fout<<"Generation";
  for (auto const& nm : momNames) fout<<","<<nm;
  fout<<"\n";

  auto engVar = engine_->toEigenVectorVariant();
  visit_same_scalar(engVar, y, [&](auto const& engVec, auto& vec){
    using Scalar = typename std::decay_t<decltype(vec)>::Scalar;
    for (size_t i=0;i<steps;++i) {
      long long gen = static_cast<long long>(startGen_) - static_cast<long long>(i)*interval;
      fout<<gen;
      for (auto idx : indices)
        fout<<","<<static_cast<double>(vec(Eigen::Index(idx)));
      fout<<"\n";
      if (i+1<steps)
        for (size_t k=0;k<interval;++k)
          vec = engVec * vec;
    }
  });

  fout.reset();
  raw.close();
}

void Epoch::testSteadyState()
{
  // optional per-operator tests
}
