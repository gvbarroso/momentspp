/*
 * Authors: Gustavo V. Barroso
 * Created: 31/08/2022
 * Last modified: 18/09/2025
 *
 */

#include <ios>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <fstream>

#include "Migration.hpp"
#include "Epoch.hpp"

void Epoch::fireParameterChanged(const bpp::ParameterList& params)
{
  if(!matchParametersValues(params))
    return;

  updateOperators_(params);

  if(operators_.empty())
    throw bpp::Exception("Epoch::attempted to fireParametersChanged in empty Operators!");

  auto acc = operators_.front()->getTransitionMatrixVariantEigen();

  for(size_t i = 1; i < operators_.size(); ++i)
  {
    auto next = operators_[i]->getTransitionMatrixVariantEigen();

    // adds next into acc in-place. Both variants must hold Eigen::SparseMatrix<Scalar>.
    std::visit([&](auto& A, const auto& B) -> void
    {
      using MatA = std::decay_t<decltype(A)>;
      using MatB = std::decay_t<decltype(B)>;
      static_assert(std::is_same_v<typename MatA::Scalar, typename MatB::Scalar>, "Matrix scalar types must match");
      A += B;
    }, acc, next);
  }

  // Hand the accumulated Eigen-variant to the engine
  engine_->setMatrix(acc);
}

void Epoch::computeExpectedSumStatsDiscrete(const MatrixEngine::VectorVariantEigen& y)
{
  auto Mvar = engine_->toEigenMatrixVariant();

  visit_eigen(y, Mvar, [&](const auto& src, auto& dst) -> void
  {
    using VecT   = std::decay_t<decltype(src)>;
    using Scalar = typename VecT::Scalar;

    if(src.size() != dst.size())
      throw bpp::Exception("computeExpectedSumStatsDiscrete: size mismatch");

    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> tmp = src;
    for(size_t k = 0; k < duration(); ++k)
      tmp = src * tmp;

    dst = tmp;
  });

  updateMoments(engine_->toEigenVectorVariant());
}

std::vector<size_t> Epoch::fetchSelectedPopIds()
{
  std::vector<size_t> ret(0);
  ret.reserve(pops_.size());

  for(auto it = std::begin(pops_); it != std::end(pops_); ++it)
  {
    if((*it)->hasSelection())
      ret.emplace_back((*it)->getId());
  }

  return ret;
}

// copies moment expectations from previous epoch into y, respecting ancestry relationships
void Epoch::transferStatistics(MatrixEngine::VectorVariantEigen& y) const
{
  auto Vvar = engine_->toEigenVectorVariant();

  visit_eigen(y, Vvar, [&](const auto& src, auto& dst) -> void
  {
    using VecT   = std::decay_t<decltype(src)>;
    using Scalar = typename VecT::Scalar;

    if(src.size() != dst.size())
      throw bpp::Exception("transferStatistics: size mismatch");

    for(Eigen::Index i = 0; i < src.size(); ++i)
      dst(i) = src(i);
  });
}

void Epoch::updateMoments(const MatrixEngine::VectorVariantEigen& y)
{
  std::visit([&](const auto& vec) -> void
  {
    using VecType = std::decay_t<decltype(vec)>;
    using Scalar  = typename VecType::Scalar;
    constexpr bool is_eigen_vector = std::is_same_v<VecType, Eigen::Matrix<Scalar, Eigen::Dynamic, 1>>;

    static_assert(is_eigen_vector, "Expected an Eigen dynamic column vector in VectorVariantEigen");

    const size_t basisSize = ssl_.getBasis().size();
    const Eigen::Index vecSize = vec.size();

    if(vecSize != static_cast<Eigen::Index>(basisSize))
      throw std::runtime_error("Epoch::updateMoments: vector size does not match ssl basis size");

    for(size_t i = 0; i < basisSize; ++i)
    {
      const Eigen::Index idx = static_cast<Eigen::Index>(i);
      double value = static_cast<double>(vec(idx));
      ssl_.getBasis()[i]->setValue(value);
    }
  }, y);
}


void Epoch::printMoments(std::ostream& stream)
{
  std::vector<std::shared_ptr<Moment>> tmp = getSslib().getBasis();

  for(auto& m : tmp)
    stream << std::setprecision(16) << m->getName() << " = " << m->getValue() << "\n";
}

void Epoch::printMomentsIntermediate(MatrixEngine::VectorVariantEigen& y,
                                     const std::string& modelName,
                                     size_t interval,
                                     const std::vector<std::string>& momNames)
{
  transferStatistics(y);

  std::string fileName = modelName + "_" + name_ + "_moments.csv.gz";
  std::ofstream rawFile(fileName, std::ios::binary);
  boost::iostreams::filtering_ostream fout;
  fout.push(boost::iostreams::gzip_compressor());
  fout.push(rawFile);

  const auto& basis = ssl_.getBasis();
  const size_t steps = duration() / interval + 1;
  std::vector<size_t> indices;
  for(size_t i = 0; i < basis.size(); ++i)
    if(std::find(momNames.begin(), momNames.end(), basis[i]->getName())
       != momNames.end())
      indices.push_back(i);

  fout << "Generation";

  for(auto& nm : momNames)
    fout << "," << nm;

  fout << "\n";

  auto VecVar = engine_->toEigenVectorVariant();
  visit_eigen(y, VecVar, [&](auto& vec, auto& engVec) -> void
  {
    using VecT   = std::decay_t<decltype(vec)>;
    for(size_t i = 0; i < steps; ++i)
    {
      fout << (startGen_ - i*interval);

      for(auto idx : indices)
        fout << "," << double(vec(Eigen::Index(idx)));

      fout << "\n";

      if(i+1 < steps)
      {
        for(size_t k = 0; k < interval; ++k)
          vec = engVec * vec;
      }
    }
  });

  fout.reset();
  rawFile.close();
}


void Epoch::printRecursions(std::ostream& stream)
{
  stream << "\n";
  const auto& basisVec = ssl_.getBasis();

  for(size_t i = 0; i < basisVec.size(); ++i)
  {
    const auto& basis = basisVec[i];
    if(basis->getName() == "I") continue;

    Eigen::Index pos = static_cast<Eigen::Index>(basis->getPosition());
    stream << "\u0394[" << basis->getName() << "] = ";

    for(auto& op : operators_)
    {
      for(size_t k = 0; k < op->getParameters().size(); ++k)
      {
        auto& param = op->getParameters()[k];
        auto matVar = op->getMatrix(k).getMatrixVariant();

        visit_eigen(matVar, [&](const auto& M) -> void
        {
          using MatT   = std::decay_t<decltype(M)>;
          using Scalar = typename MatT::Scalar;
          Scalar scale = Scalar(param.getValue());
          if(scale == Scalar(0)) return;

          for(Eigen::Index l = 0; l < M.cols(); ++l)
          {
            Scalar c = M.coeff(pos,l) / scale;

            if(c == Scalar(0))
              continue;

            if(c>Scalar(0))
              stream << "+";

            stream << std::fixed << std::setprecision(3) << c << "*" << param.getName() << "*" << ssl_.getBasis()[size_t(l)]->getName() << " ";
          }
        });
      }
    }
    stream << "\n";
  }
}


void Epoch::printTransitionMat(const std::string& fileName) const
{
  auto Mvar = engine_->toEigenMatrixVariant();
  visit_eigen(Mvar, [&](const auto& M) -> void
  {
    std::ofstream out(fileName, std::ios::binary);

    if(!out)
      throw std::runtime_error("Cannot open " + fileName);

    out << std::scientific << std::setprecision(12) << M.rows() << " " << M.cols() << " " << M.nonZeros() << "\n";

    for(Eigen::Index k = 0; k < M.outerSize(); ++k)
    {
      for(typename std::decay_t<decltype(M)>::InnerIterator it(M, k); it; ++it)
        out << it.row() << " " << it.col() << " " << static_cast<long double>(it.value()) << "\n";
    }
  });
}


void Epoch::computeEigenSteadyState()
{
  EigenResult res = findLeadingEigenpair();

  if(res.value.toDouble() > 1.0 + 1e-14)
  {
    std::cout << "Condition Number = " << fetchConditionNumber() << "\n" << "Leading eigenvalue = " << res.value.toDouble() << "\n";
    throw bpp::Exception("Leading eigenvalue > 1");
  }

  auto vecVar = engine_->toEigenVectorVariant();
  visit_eigen(vecVar, [&](auto& dst) -> void
  {
    using VecT   = std::decay_t<decltype(dst)>;
    using Scalar = typename VecT::Scalar;

    if(dst.size() != static_cast<Eigen::Index>(res.vector.size()))
      throw bpp::Exception("computeEigenSteadyState: size mismatch");

    for(Eigen::Index i = 0; i < dst.size(); ++i)
      dst(i) = static_cast<Scalar>(res.vector(i));
  });
  engine_->setVectorFromEigen(std::move(vecVar));

  updateMoments(engine_->toEigenVectorVariant());
}

// computePseudoSteadyStateDiscrete
void Epoch::computePseudoSteadyStateDiscrete(double tol)
{
  bool converged = false;
  auto Mvar = engine_->toEigenMatrixVariant();

  visit_eigen(Mvar, [&](const auto& A) -> void
  {
    using Scalar = typename std::decay_t<decltype(A)>::Scalar;
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> y(A.cols());

    // parameters for initialization
    size_t pop   = ssl_.getPopIndices()[0];
    double mu    = getParameterValue("u_"  + bpp::TextTools::toString(pop));
    double s     = getParameterValue("s_"  + bpp::TextTools::toString(pop));
    size_t twoN  = static_cast<size_t>(1.0 / getParameterValue("1/2N_" + bpp::TextTools::toString(pop)));

    double hr = twoN * mu;
    double hl = 2.0 * twoN * mu * std::exp(2.0 * twoN * s) / (std::exp(2.0 * twoN * s) - 1.0) - mu / s;
    double f  = 1.0;

    // initial guess
    for(Eigen::Index i = 0; i < y.size(); ++i)
    {
      const std::string& prefix = ssl_.getBasis()[static_cast<size_t>(i)]->getPrefix();
      if(prefix == "Hl")       y(i) = Scalar(hl * f);
      else if(prefix == "Hr")  y(i) = Scalar(hr * f);
      else if(prefix == "pi2") y(i) = Scalar(hr * hl * f);
      else if(prefix == "I")   y(i) = Scalar(1.0);
      else                     y(i) = Scalar(hr * hl * f * 1e-1);

      if(i > 0 && prefix == ssl_.getBasis()[static_cast<size_t>(i - 1)]->getPrefix())
        f *= 0.925;

      else
        f = 1.0;
    }

    // burn-in iterations
    for(size_t b = 0; b < twoN / 10; ++b)
      y = A * y;

    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> prev = y;
    size_t maxIter = 20 * twoN;

    // power-method iterations until tolerance
    for(size_t iter = 0; iter < maxIter; ++iter)
    {
      y = A * y;

      double maxRel = 0.0;
      for(Eigen::Index i = 0; i < y.size(); ++i)
      {
        double pv = static_cast<double>(prev(i));
        double cv = static_cast<double>(y(i));
        double rel = std::abs(cv - pv) / std::max(1.0, std::abs(pv));
        maxRel = std::max(maxRel, rel);
      }

      if(maxRel <= tol)
      {
        std::cout << "Pseudo steady-state converged after "
                  << iter << " iterations, " << name_ << "\n";
        converged = true;
        break;
      }

      prev = y;
    }

    if(converged)
    {
      auto vecVar = engine_->toEigenVectorVariant();
      visit_eigen(vecVar, [&](auto& dst) -> void { dst = y; });
      engine_->setVectorFromEigen(std::move(vecVar));
    }
  });

  if(!converged)
  {
    std::cerr << "Pseudo steady-state did not converge. "
              << "Falling back to eigen steady-state.\n";
    computeEigenSteadyState();
    return;
  }

  updateMoments(engine_->toEigenVectorVariant());
}

// computePseudoSteadyStateContinuous
void Epoch::computePseudoSteadyStateContinuous(double burnInTime,
                                               double dt,
                                               double tol)
{
  bool converged = false;
  auto Mvar = engine_->toEigenMatrixVariant();

  visit_eigen(Mvar, [&](const auto& A) -> void
  {
    using Scalar = typename std::decay_t<decltype(A)>::Scalar;
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> y(A.cols());

    // parameters for initialization
    size_t pop   = ssl_.getPopIndices()[0];
    double mu    = getParameterValue("u_"  + bpp::TextTools::toString(pop));
    double s     = getParameterValue("s_"  + bpp::TextTools::toString(pop));
    size_t twoN  = static_cast<size_t>(1.0 / getParameterValue("1/2N_" + bpp::TextTools::toString(pop)));

    double hr = twoN * mu;
    double hl = 2.0 * twoN * mu * std::exp(2.0 * twoN * s) / (std::exp(2.0 * twoN * s) - 1.0) - mu / s;
    double f  = 1.0;

    // initial guess
    for(Eigen::Index i = 0; i < y.size(); ++i)
    {
      const std::string& prefix = ssl_.getBasis()[static_cast<size_t>(i)]->getPrefix();
      if(prefix == "Hl")       y(i) = Scalar(hl * f);
      else if(prefix == "Hr")  y(i) = Scalar(hr * f);
      else if(prefix == "pi2") y(i) = Scalar(hr * hl * f);
      else if(prefix == "I")   y(i) = Scalar(1.0);
      else                     y(i) = Scalar(hr * hl * f * 1e-1);

      if(i > 0 && prefix == ssl_.getBasis()[static_cast<size_t>(i - 1)]->getPrefix())
        f *= 0.925;

      else
        f = 1.0;
    }

    // burn-in via integrateTyped single-step evolution
    size_t burnSteps = static_cast<size_t>(burnInTime / dt);
    for(size_t b = 0; b < burnSteps; ++b)
      y = integrateTyped<Scalar>(y, dt, dt);

    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> prev = y;
    size_t maxSteps = static_cast<size_t>(10 * burnInTime / dt);

    // iterative evolution until convergence
    for(size_t iter = 0; iter < maxSteps; ++iter)
    {
      prev = y;
      y = integrateTyped<Scalar>(y, dt, dt);

      double maxRel = 0.0;

      for(Eigen::Index i = 0; i < y.size(); ++i)
      {
        double pv = static_cast<double>(prev(i));
        double cv = static_cast<double>(y(i));
        double rel = std::abs(cv - pv) / std::max(1.0, std::abs(pv));
        maxRel = std::max(maxRel, rel);
      }

      if(maxRel <= tol)
      {
        std::cout << "Continuous pseudo steady-state converged after "
                  << iter << " steps, " << name_ << "\n";
        converged = true;
        break;
      }
    }

    if(converged)
    {
      auto vecVar = engine_->toEigenVectorVariant();
      visit_eigen(vecVar, [&](auto& dst) -> void { dst = y; });
      engine_->setVectorFromEigen(std::move(vecVar));
    }
  });

  if(!converged)
  {
    std::cerr << "Crank-Nicolson steady-state did not converge. "
              << "Falling back to eigen-based steady-state.\n";
    computeEigenSteadyState();
    return;
  }

  updateMoments(engine_->toEigenVectorVariant());
}

// test existence of steady-state in models with gene-flow
void Epoch::testSteadyState()
{
  /*
  if(pops_.size() > 1)
  {
    for(size_t i = 0; i < operators_.size(); ++i)
    {
      auto tmp = std::dynamic_pointer_cast<Migration>(operators_[i]);

      if(tmp != nullptr)
        tmp->testFlow();
    }
  }
  */
}

template <typename Scalar>
Eigen::Matrix<Scalar, Eigen::Dynamic, 1>
Epoch::integrateTyped(const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& moms, double dt, double totalTime) const
{
  const Eigen::Index steps = static_cast<Eigen::Index>(
    std::max<int>(1, static_cast<int>(std::ceil(totalTime / dt)))
  );
  Eigen::Matrix<Scalar, Eigen::Dynamic, 1> y = moms;

  auto Avar = engine_->toEigenMatrixVariant();
  visit_eigen(Avar, [&](const auto& A) -> void
  {
    for(Eigen::Index t = 0; t < steps; ++t)
      y = A * y;
  });

  return y;
}

MatrixEngine::VectorVariantEigen Epoch::integrate(const MatrixEngine::VectorVariantEigen& moms, double dt, double totalTime) const
{
  return std::visit([&](const auto& vec) -> MatrixEngine::VectorVariantEigen
  {
    using VecT = std::decay_t<decltype(vec)>;
    using Scalar = typename VecT::Scalar;

    auto res = integrateTyped<Scalar>(vec, dt, totalTime);
    return MatrixEngine::VectorVariantEigen(std::move(res));
  }, moms);
}

template <typename Scalar>
Eigen::Matrix<Scalar, Eigen::Dynamic, 1>
Epoch::integrateAdaptiveTyped(const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& moms,
                              double dt,
                              double totalTime,
                              double tolerance,
                              double dtMin,
                              double dtMax) const
{
  auto Avar = engine_->toEigenMatrixVariant();
  Eigen::Matrix<Scalar, Eigen::Dynamic, 1> y = moms;
  double t = 0.0;
  double h = dt;

  visit_eigen(Avar, [&](const auto& A) -> void
  {
    while(t < totalTime)
    {
      h = std::clamp(h, dtMin, dtMax);
      if(t + h > totalTime) h = totalTime - t;

      Eigen::Matrix<Scalar, Eigen::Dynamic, 1> y1 = A * y;
      Eigen::Matrix<Scalar, Eigen::Dynamic, 1> y2 = A * (A * y);

      double err = 0.0;
      for(Eigen::Index i = 0; i < y.size(); ++i)
        err = std::max(err, std::abs(double(y1(i) - y2(i))));

      if(err <= tolerance)
      {
        y = y2;
        t += h;
        h = std::min(h * 1.5, dtMax);
      }

      else
        h = std::max(h * 0.5, dtMin);
    }
  });

  return y;
}

MatrixEngine::VectorVariantEigen Epoch::integrateAdaptive(const MatrixEngine::VectorVariantEigen& moms, double dt, double totalTime, double tolerance) const
{
  return std::visit([&](const auto& vec) -> MatrixEngine::VectorVariantEigen
  {
    using VecT   = std::decay_t<decltype(vec)>;
    using Scalar = typename VecT::Scalar;

    auto res = integrateAdaptiveTyped<Scalar>(vec, dt, totalTime, tolerance);
    return MatrixEngine::VectorVariantEigen(std::move(res));
  }, moms);
}

void Epoch::init_()
{
  if(operators_.empty())
    throw bpp::Exception("Epoch::init_() called with no operators.");

  if(!engine_)
    throw bpp::Exception("Epoch::init_() called with null engine_.");

  // Start with the first operator's Eigen-based transition variant
  auto acc = operators_.front()->getTransitionMatrixVariantEigen();

  for(size_t i = 1; i < operators_.size(); ++i) {
    auto next = operators_[i]->getTransitionMatrixVariantEigen();

    std::visit([&](auto& A, const auto& B) -> void
    {
      using MatA = std::decay_t<decltype(A)>;
      using MatB = std::decay_t<decltype(B)>;
      static_assert(std::is_same_v<typename MatA::Scalar, typename MatB::Scalar>,
                    "Matrix scalar types must match when accumulating transitions");
      A += B;
    }, acc, next);
  }

  // Commit accumulated Eigen variant to the engine in one call
  engine_->setMatrix(acc);

  // turn "delta" into full transition, then tidy up
  engine_->addIdentityInPlace();
  engine_->pruneInPlace();
  engine_->compressInPlace();

  // testSteadyState();
}
