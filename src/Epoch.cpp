/*
 * Authors: Gustavo V. Barroso
 * Created: 31/08/2022
 * Last modified: 20/09/2025
 *
 */

#include <ios>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <fstream>
#include <variant>
#include <type_traits>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <fstream>

#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/SparseLU>  // for double‐type LU
#include <eigen3/Eigen/LU>

#include "Migration.hpp"
#include "Epoch.hpp"

void Epoch::fireParameterChanged(const bpp::ParameterList& params)
{
  if (!matchParametersValues(params))
    return;

  for (auto& op : operators_)
    op->fireParameterChanged(params);

  init_();
}

void Epoch::computeExpectedSumStatsDiscrete(const MatrixEngine::VectorVariantEigen& y)
{
  auto matVar = engine_->toEigenMatrixVariant();
  MatrixEngine::VectorVariantEigen result;

  visit_same_scalar(matVar, y, [&](auto const& A, auto const& vec) {
    using Scalar = typename std::decay_t<decltype(vec)>::Scalar;

    if (vec.size() != A.cols())
      throw bpp::Exception("computeExpectedSumStatsDiscrete: size mismatch");

    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> tmp = vec;
    for (size_t k = 0; k < duration(); ++k)
      tmp = A * tmp;

    result = std::move(tmp);
  });

  updateMoments(result);
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
// Transfer a steady-state from a previous Epoch’s vector into this one.
// prevY must be a VectorVariantEigen holding either a VectorXd or an mpreal vector.
void Epoch::transferStatistics(const MatrixEngine::VectorVariantEigen& prevY)
{
  // 1) Allocate a new Eigen‐vector of the correct type & size
  auto newY = engine_->toEigenVectorVariant();
  const size_t n = ssl_.getBasis().size();

  // 2) Remap: for each i in [0..n-1], pick the parent position from prevY
  std::visit([&](auto& dst, auto const& src) {
    using DVec    = std::decay_t<decltype(dst)>;
    using SScalar = typename std::decay_t<decltype(src)>::Scalar;
    using DScalar = typename DVec::Scalar;

    // Ensure we’re not mixing doubles with mpreals
    static_assert(std::is_same_v<SScalar, DScalar>,
                  "Epoch::transferStatistics: scalar types must match");

    if (static_cast<size_t>(dst.size()) != n)
      throw bpp::Exception("Epoch::transferStatistics: target size mismatch");

    for (size_t i = 0; i < n; ++i)
    {
      // Find the parent‐index for basis[i]
      size_t parentPos =
        ssl_.getBasis()[i]->getParent()->getPosition();

      if (parentPos >= static_cast<size_t>(src.size()))
        throw bpp::Exception("Epoch::transferStatistics: parent index out of range");

      // Copy that entry
      dst(Eigen::Index(i)) = src(Eigen::Index(parentPos));
    }
  }, newY, prevY);

  // 3) Store the remapped vector back into the engine
  engine_->setVectorFromEigen(std::move(newY));
}

void Epoch::updateMoments(const MatrixEngine::VectorVariantEigen& y)
{
  std::visit([&](auto const& vec) {
    for (size_t i = 0; i < ssl_.getBasis().size(); ++i)
      ssl_.getBasis()[i]->setValue(static_cast<double>(vec(Eigen::Index(i))));
  }, y);
}

void Epoch::printMoments(std::ostream& stream)
{
  std::vector<std::shared_ptr<Moment>> tmp = getSslib().getBasis();

  for(auto& m : tmp)
    stream << std::setprecision(16) << m->getName() << " = " << m->getValue() << "\n";
}

//------------------------------------------------------------------------------
// Write out CSV.gz of selected moments over time, stepping by 'interval'
// yPrev is a VectorVariantEigen from the previous epoch.
//------------------------------------------------------------------------------
void Epoch::printMomentsIntermediate(
    MatrixEngine::VectorVariantEigen& yPrev,
    const std::string& modelName,
    size_t interval,
    const std::vector<std::string>& momNames)
{
    // 1) Remap the incoming vector into this epoch’s state
    transferStatistics(yPrev);

    // 2) Prepare compressed output file
    std::string fileName = modelName + "_" + name_ + "_moments.csv.gz";
    std::ofstream rawFile(fileName, std::ios::binary);
    if (!rawFile)
        throw bpp::Exception("printMomentsIntermediate: cannot open " + fileName);

    boost::iostreams::filtering_ostream fout;
    fout.push(boost::iostreams::gzip_compressor());
    fout.push(rawFile);

    // 3) Decide which basis indices to record
    const auto& basis = ssl_.getBasis();
    size_t nStates = basis.size();
    size_t steps   = duration() / interval + 1;

    std::vector<size_t> indices;
    indices.reserve(momNames.size());
    for (size_t i = 0; i < nStates; ++i)
        if (std::find(momNames.begin(), momNames.end(), basis[i]->getName())
            != momNames.end())
            indices.push_back(i);

    // 4) Write CSV header
    fout << "Generation";
    for (auto const& nm : momNames) fout << "," << nm;
    fout << "\n";

    // 5) Grab the transition matrix and the remapped initial vector
    auto matVar = engine_->toEigenMatrixVariant();
    auto vecVar = engine_->toEigenVectorVariant();

    // 6) Iterate, update by applying M^interval, and dump each snapshot
    visit_same_scalar(matVar, vecVar, [&](auto const& M, auto& v) {
        using Scalar = typename std::decay_t<decltype(v)>::Scalar;

        for (size_t step = 0; step < steps; ++step)
        {
            long long generation =
              static_cast<long long>(startGen_)
              - static_cast<long long>(step) * static_cast<long long>(interval);

            // Print this generation’s values
            fout << generation;
            for (size_t idx : indices)
                fout << "," << static_cast<double>(v(Eigen::Index(idx)));
            fout << "\n";

            // Advance by 'interval' applications of M
            if (step + 1 < steps)
            {
                for (size_t k = 0; k < interval; ++k)
                    v = M * v;
            }
        }
    });

    // 7) Clean up
    fout.reset();
    rawFile.close();
}

//------------------------------------------------------------------------------
// Print recursion formulas Δ[state] = … for each basis function except “I”
//------------------------------------------------------------------------------
void Epoch::printRecursions(std::ostream& stream)
{
    stream << '\n';
    const auto& basis = ssl_.getBasis();

    for (size_t i = 0; i < basis.size(); ++i)
    {
        auto const& state = basis[i];
        if (state->getName() == "I")
            continue;

        Eigen::Index pos = state->getPosition();
        stream << "\u0394[" << state->getName() << "] = ";

        bool firstTerm = true;

        for (auto const& op : operators_)
        {
            auto const& params = op->getParameters();

            for (size_t k = 0; k < params.size(); ++k)
            {
                auto const& param = params[k];
                double pval = param.getValue();
                if (pval == 0.0)
                    continue;

                // 1) get the operator’s Matrix<Scalar> variant
                auto matVar = op->getMatrix(k).getMatrixVariant();

                // 2) convert to our Eigen‐variant
                MatrixEngine::MatrixVariantEigen eigenVar =
                  std::visit([](auto const& Mw) {
                    return MatrixEngine::MatrixVariantEigen{ Mw.eigen() };
                  }, matVar);

                // 3) visit by scalar type and print nonzero terms
                std::visit([&](auto const& M) {
                    using MatT    = std::decay_t<decltype(M)>;
                    using Scalar  = typename MatT::Scalar;

                    Scalar scale = static_cast<Scalar>(pval);
                    Eigen::Index ncols = M.cols();

                    for (Eigen::Index col = 0; col < ncols; ++col)
                    {
                        Scalar coeff = M.coeff(pos, col) / scale;
                        if (coeff == Scalar(0))
                            continue;

                        if (!firstTerm && coeff > Scalar(0))
                            stream << '+';

                        stream << std::fixed << std::setprecision(3)
                               << coeff
                               << '*' << param.getName()
                               << '*' << basis[size_t(col)]->getName()
                               << ' ';

                        firstTerm = false;
                    }
                }, eigenVar);
            }
        }

        stream << '\n';
    }
}

//------------------------------------------------------------------------------
// Dump the (summed) transition matrix to a file as row col value triples
//------------------------------------------------------------------------------
void Epoch::printTransitionMat(const std::string& fileName) const
{
    auto eigenVar = engine_->toEigenMatrixVariant();

    std::ofstream out(fileName, std::ios::binary);
    if (!out)
        throw bpp::Exception("printTransitionMat: cannot open " + fileName);

    // preserve caller’s formatting
    auto oldFlags = out.flags();
    auto oldPrec  = out.precision();

    out << std::scientific << std::setprecision(16);

    std::visit([&](auto const& M) {
        using MatT   = std::decay_t<decltype(M)>;
        using Scalar = typename MatT::Scalar;

        out << M.rows() << ' ' << M.cols() << ' ' << M.nonZeros() << '\n';

        for (Eigen::Index k = 0; k < M.outerSize(); ++k)
        {
            for (typename MatT::InnerIterator it(M, k); it; ++it)
            {
                out << it.row() << ' '
                    << it.col() << ' '
                    << it.value() << '\n';
            }
        }
    }, eigenVar);

    // restore formatting
    out.flags(oldFlags);
    out.precision(oldPrec);
}

void Epoch::computeEigenSteadyState()
{
  EigenResult res = findLeadingEigenpair();

  auto matVar = engine_->toEigenMatrixVariant();
  auto vecVar = std::visit([&](auto const& M) {
    using Scalar = typename std::decay_t<decltype(M)>::Scalar;
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> v(M.cols());
    for (Eigen::Index i = 0; i < M.cols(); ++i)
      v(i) = static_cast<Scalar>(res.vector(i));
    return MatrixEngine::VectorVariantEigen(std::move(v));
  }, matVar);

  engine_->setVectorFromEigen(std::move(vecVar));
  updateMoments(engine_->toEigenVectorVariant());
}

void Epoch::computePseudoSteadyStateDiscrete(double tol)
{
  bool converged = false;
  auto matVar = engine_->toEigenMatrixVariant();

  // 1) Power‐method per scalar type, with custom init & burn‐in
  std::visit([&](auto const& A) {
    using MatT    = std::decay_t<decltype(A)>;
    using Scalar  = typename MatT::Scalar;
    using Vec     = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

    // --- build and initialize y ---
    Vec y(A.cols());

    // grab pop, μ, s, 2N from your SumStatsLibrary & params
    size_t pop   = ssl_.getPopIndices()[0];
    double mu    = getParameterValue("u_"   + bpp::TextTools::toString(pop));
    double s     = getParameterValue("s_"   + bpp::TextTools::toString(pop));
    size_t twoN  = static_cast<size_t>(1.0 /
                   getParameterValue("1/2N_" + bpp::TextTools::toString(pop)));

    // compute hr, hl
    double hr = twoN * mu;
    double hl;
    if (std::abs(s) < 1e-14)
      hl = twoN * mu;
    else
      hl = 2.0 * twoN * mu * std::exp(2.0*twoN*s)
           / (std::exp(2.0*twoN*s) - 1.0)
         - mu / s;

    double f = 1.0;
    for (Eigen::Index i = 0; i < y.size(); ++i)
    {
      const auto& prefix = ssl_.getBasis()[size_t(i)]->getPrefix();
      if      (prefix == "Hl")   y(i) = Scalar(hl * f);
      else if (prefix == "Hr")   y(i) = Scalar(hr * f);
      else if (prefix == "pi2")  y(i) = Scalar(hr * hl * f);
      else if (prefix == "I")    y(i) = Scalar(1.0);
      else                       y(i) = Scalar(hr * hl * f * 1e-1);

      if (i > 0 &&
          prefix == ssl_.getBasis()[size_t(i - 1)]->getPrefix())
        f *= 0.925;
      else
        f = 1.0;
    }

    // --- burn‐in iterations ---
    size_t burnSteps = twoN / 10;
    for (size_t b = 0; b < burnSteps; ++b)
      y = A * y;

    // --- power‐method until convergence ---
    Vec prev = y;
    size_t maxIter = 20 * twoN;
    for (size_t iter = 0; iter < maxIter; ++iter)
    {
      y = A * y;
      double maxRel = 0.0;
      for (Eigen::Index i = 0; i < y.size(); ++i)
      {
        double pv = static_cast<double>(prev(i));
        double cv = static_cast<double>(y(i));
        double rel = std::abs(cv - pv) / std::max(1.0, std::abs(pv));
        maxRel = std::max(maxRel, rel);
      }
      if (maxRel <= tol)
      {
        std::cout
          << "Pseudo steady-state converged after "
          << iter << " iterations, " << name_ << "\n";
        converged = true;
        break;
      }
      prev = y;
    }

    // 2) commit if converged
    if (converged)
    {
      auto vecVar = engine_->toEigenVectorVariant();
      std::visit([&](auto& dst) {
        using VecT    = std::decay_t<decltype(dst)>;
        using DScalar = typename VecT::Scalar;
        // compile‐time guard
        if constexpr (!std::is_same_v<DScalar, Scalar>)
          throw bpp::Exception("Type mismatch in pseudo‐steady assignment");
        dst = y;
      }, vecVar);
      engine_->setVectorFromEigen(std::move(vecVar));
    }
  }, matVar);

  // 3) fallback on non‐convergence
  if (!converged)
  {
    std::cerr
      << "Pseudo steady-state did not converge. "
      << "Falling back to eigen steady-state.\n";
    computeEigenSteadyState();
    return;
  }

  // 4) update moments
  updateMoments(engine_->toEigenVectorVariant());
}

//------------------------------------------------------------------------------
// Typed helper: does all the work in the Scalar domain
// Returns true if converged, and writes the steady-state vector into outY
//------------------------------------------------------------------------------
template<typename Scalar>
bool Epoch::computePseudoSSContinuousTyped(
    const Eigen::SparseMatrix<Scalar>& A,
    double burnInTime,
    double dt,
    double tol,
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& outY) const
{
  using Vec   = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
  using MWrap = Matrix<Scalar>;
  using VWrap = Vector<Scalar>;

  const Eigen::Index n = A.cols();
  Vec y(n);

  // 1) Rough initial guess (hr, hl, prefix‐based)
  size_t pop   = ssl_.getPopIndices()[0];
  double mu    = getParameterValue("u_"   + bpp::TextTools::toString(pop));
  double s     = getParameterValue("s_"   + bpp::TextTools::toString(pop));
  size_t twoN  = static_cast<size_t>(1.0 /
                   getParameterValue("1/2N_" + bpp::TextTools::toString(pop)));

  double hr = twoN * mu;
  double hl;
  if (std::abs(s) < 1e-14)
    hl = twoN * mu;
  else
    hl = 2.0*twoN*mu*std::exp(2.0*twoN*s)
         / (std::exp(2.0*twoN*s) - 1.0)
       - mu/s;

  double f = 1.0;
  for (Eigen::Index i = 0; i < n; ++i)
  {
    const auto& prefix = ssl_.getBasis()[size_t(i)]->getPrefix();
    if      (prefix == "Hl")   y(i) = Scalar(hl * f);
    else if (prefix == "Hr")   y(i) = Scalar(hr * f);
    else if (prefix == "pi2")  y(i) = Scalar(hr * hl * f);
    else if (prefix == "I")    y(i) = Scalar(1.0);
    else                       y(i) = Scalar(hr * hl * f * 1e-1);

    if (i > 0 &&
        prefix == ssl_.getBasis()[size_t(i - 1)]->getPrefix())
      f *= 0.925;
    else
      f = 1.0;
  }

  // 2) Burn-in: repeated single-step CN
  size_t burnSteps = static_cast<size_t>(burnInTime / dt);
  MWrap Mwrap(A);
  VWrap Vwrap;
  for (size_t b = 0; b < burnSteps; ++b)
  {
    Vwrap.eigen() = y;
    if constexpr (std::is_same_v<Scalar, double>)
      y = integrateDoubleCN(Mwrap, Vwrap, dt, dt).eigen();
    else
      y = integrateMpfrCN  (Mwrap, Vwrap, dt, dt).eigen();
  }

  // 3) Fixed-step power-iteration CN until convergence
  Vec prev = y;
  size_t maxSteps = static_cast<size_t>(10 * burnInTime / dt);
  bool converged = false;

  for (size_t iter = 0; iter < maxSteps; ++iter)
  {
    Vwrap.eigen() = y;
    if constexpr (std::is_same_v<Scalar, double>)
      y = integrateDoubleCN(Mwrap, Vwrap, dt, dt).eigen();
    else
      y = integrateMpfrCN  (Mwrap, Vwrap, dt, dt).eigen();

    double maxRel = 0.0;
    for (Eigen::Index i = 0; i < n; ++i)
    {
      double pv = static_cast<double>(prev(i));
      double cv = static_cast<double>(y(i));
      double rel = std::abs(cv - pv) / std::max(1.0, std::abs(pv));
      maxRel = std::max(maxRel, rel);
    }

    if (maxRel <= tol)
    {
      std::cout
        << "Continuous pseudo steady-state converged after "
        << iter << " steps, " << name_ << "\n";
      converged = true;
      break;
    }
    prev = y;
  }

  if (!converged)
    return false;

  // 4) Commit outY
  outY = std::move(y);
  return true;
}

//------------------------------------------------------------------------------
// Public dispatcher: unpack the variant and call the typed helper
//------------------------------------------------------------------------------
void Epoch::computePseudoSteadyStateContinuous(
    double burnInTime,
    double dt,
    double tol)
{
  auto matVar = engine_->toEigenMatrixVariant();
  bool converged = false;

  // Prepare storage for each scalar
  Eigen::Matrix<double,      Eigen::Dynamic, 1> yD;
  Eigen::Matrix<mpfr::mpreal, Eigen::Dynamic, 1> yM;

  std::visit([&](auto const& A) {
    using Scalar = typename std::decay_t<decltype(A)>::Scalar;

    if constexpr (std::is_same_v<Scalar, double>)
      converged = computePseudoSSContinuousTyped(A, burnInTime, dt, tol, yD);
    else
      converged = computePseudoSSContinuousTyped(A, burnInTime, dt, tol, yM);
  }, matVar);

  if (!converged)
  {
    std::cerr
      << "Crank–Nicolson steady-state did not converge. "
      << "Falling back to eigen-based steady-state.\n";
    computeEigenSteadyState();
    return;
  }

  // 5) Build & commit the vector variant
  MatrixEngine::VectorVariantEigen vecVar =
    (yD.size() ? VectorVariantEigen(std::move(yD))
               : VectorVariantEigen(std::move(yM)));

  engine_->setVectorFromEigen(std::move(vecVar));
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

MatrixEngine::VectorVariant
Epoch::integrate(double dt, double totalTime) const
{
    if (dt <= 0.0 || totalTime < 0.0)
        throw bpp::Exception("Epoch::integrate: dt must be >0 and totalTime >=0");

    const auto& matVar = engine_->getMatrixVariant();
    const auto& vecVar = engine_->getVectorVariant();
    MatrixEngine::VectorVariant result;

    visit_same_scalar(matVar, vecVar, [&](auto const& M, auto const& V) {
        using Scalar = typename std::decay_t<decltype(M)>::Scalar;

        // Cast time parameters into Scalar if needed
        Scalar dtS  = static_cast<Scalar>(dt);
        Scalar Ttot = static_cast<Scalar>(totalTime);

        if constexpr (std::is_same_v<Scalar, double>) {
            result = integrateDoubleCN(M, V, dt, totalTime);
        } else {
            result = integrateMpfrCN(M, V, dtS, Ttot);
        }
    });

    return result;
}

MatrixEngine::VectorVariant
Epoch::integrateAdaptive(double dt,
                         double totalTime,
                         double tol,
                         double dtMin,
                         double dtMax) const
{
    if (dtMin <= 0.0 || dtMax <= 0.0 || dtMin > dtMax)
        throw bpp::Exception("Epoch::integrateAdaptive: invalid dtMin/dtMax");
    if (tol < 0.0)
        throw bpp::Exception("Epoch::integrateAdaptive: tol must be >=0");

    const auto& matVar = engine_->getMatrixVariant();
    const auto& vecVar = engine_->getVectorVariant();
    MatrixEngine::VectorVariant result;

    visit_same_scalar(matVar, vecVar, [&](auto const& M, auto const& V) {
        using Scalar = typename std::decay_t<decltype(M)>::Scalar;

        // Cast time args into Scalar
        Scalar dtS    = static_cast<Scalar>(dt);
        Scalar Ttot   = static_cast<Scalar>(totalTime);
        Scalar tolS   = static_cast<Scalar>(tol);
        Scalar dtMinS = static_cast<Scalar>(dtMin);
        Scalar dtMaxS = static_cast<Scalar>(dtMax);

        if constexpr (std::is_same_v<Scalar, double>) {
            result = integrateAdaptiveDoubleCN(
                M, V, dt, totalTime, tol, dtMin, dtMax
            );
        } else {
            result = integrateAdaptiveMpfrCN(
                M, V, dtS, Ttot, tolS, dtMinS, dtMaxS
            );
        }
    });

    return result;
}

void Epoch::init_()
{
  if (operators_.empty())
    throw bpp::Exception("Epoch::init_() called with no operators.");

  if (!engine_)
    throw bpp::Exception("Epoch::init_() called with null engine_.");

  auto acc = operators_.front()->getTransitionMatrixVariantEigen();
  for (size_t i = 1; i < operators_.size(); ++i) {
    auto next = operators_[i]->getTransitionMatrixVariantEigen();
    visit_same_scalar(acc, next, [](auto& A, const auto& B) {
      A += B;
    });
  }

  engine_->setMatrixFromEigen(std::move(acc));
  engine_->addIdentityInPlace();
  engine_->pruneInPlace();
  engine_->compressInPlace();
}

