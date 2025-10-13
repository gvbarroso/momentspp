/*
 * Authors: Gustavo V. Barroso
 * Created: 31/08/2022
 * Last modified: 13/10/2025
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

#include <cxxabi.h>

#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/SparseLU>
#include <eigen3/Eigen/LU>

#include "Migration.hpp"
#include "Epoch.hpp"

#include "Epoch.hpp"
#include "VariantUtils.hpp"

//------------------------------------------------------------------------------
// Responds to parameter updates by propagating to each operator and
// re-initializing this epoch.
//------------------------------------------------------------------------------
void Epoch::fireParameterChanged(const bpp::ParameterList& params)
{
  if(!matchParametersValues(params))
    return;

  for(auto& op : operators_)
    op->fireParameterChanged(params);

  init_();
}

//------------------------------------------------------------------------------
// Computes E[Y(t)] = A^duration * y0 for discrete‐time model.
//------------------------------------------------------------------------------
void Epoch::computeExpectedSumStatsDiscrete(MatrixEngine::VectorVariant& y)
{
  auto matVar = engine_->getMatrixVariant();

  visitMatrixAndVector(matVar, y, [&](auto const& A, auto& vecWrap)
  {
    for(size_t k = 0; k < duration(); ++k)
      vecWrap = A * vecWrap;
  });

  updateMoments(y);
}

//------------------------------------------------------------------------------
// Return IDs of populations under selection in this epoch.
//------------------------------------------------------------------------------
std::vector<size_t> Epoch::fetchSelectedPopIds()
{
  std::vector<size_t> ret;
  ret.reserve(pops_.size());

  for(auto it = std::begin(pops_); it != std::end(pops_); ++it)
  {
    if((*it)->hasSelection())
      ret.emplace_back((*it)->getId());
  }

  return ret;
}

//------------------------------------------------------------------------------
// Remaps moment‐expectation vector according to population ancestry relationships.
//------------------------------------------------------------------------------
void Epoch::transferStatistics(MatrixEngine::VectorVariant& prevY)
{
  MatrixEngine::VectorVariant newY = engine_->getVectorVariant();
  size_t numStatsThisEpoch = ssl_.getBasis().size();

  visitSameType(newY, prevY, [&](auto& dstWrap, auto const& srcWrap)
  {
    auto& dst = dstWrap.eigen();
    const auto& src = srcWrap.eigen();

    // for each Moment in *this Epoch, we assign its value from its parental Moment from the previous Epoch
    for(size_t i = 0; i < numStatsThisEpoch; ++i)
    {
      size_t parentPos = ssl_.getBasis()[i]->getParent()->getPosition();

      if(parentPos >= static_cast<size_t>(src.size()))
        throw bpp::Exception("transferStatistics: out-of-range");

      dst(Eigen::Index(i)) = src(Eigen::Index(parentPos));
    }
  });

  prevY = newY;
  //engine_->setVector(std::move(newY));
}

//------------------------------------------------------------------------------
// Updates the SumStatsLibrary’s basis moments from a variant Eigen vector.
//------------------------------------------------------------------------------
void Epoch::updateMoments(const MatrixEngine::VectorVariant& y)
{
  std::visit([&](auto const& wrap)
  {
    auto const& raw = wrap.eigen();

    for(size_t i = 0; i < ssl_.getBasis().size(); ++i)
      ssl_.getBasis()[i]->setValue(static_cast<double>(raw(Eigen::Index(i))));

  }, y);
}

//------------------------------------------------------------------------------
// Print the current basis moments to a text stream.
//------------------------------------------------------------------------------
void Epoch::printMoments(std::ostream& stream)
{
  auto tmp = ssl_.getBasis();

  for(auto& m : tmp)
    stream << std::setprecision(16) << m->getName() << " = " << m->getValue() << "\n";
}

//------------------------------------------------------------------------------
// Write out CSV.gz of selected moments over time, stepping by 'interval'.
//------------------------------------------------------------------------------
void Epoch::printMomentsIntermediate(MatrixEngine::VectorVariant& yPrev,
                                     const std::string& modelName,
                                     size_t interval,
                                     const std::vector<std::string>& momNames)
{
  transferStatistics(yPrev);

  std::string fileName = modelName + "_" + name_ + "_moments.csv.gz";
  std::ofstream rawFile(fileName, std::ios::binary);

  if(!rawFile)
    throw bpp::Exception("printMomentsIntermediate: cannot open " + fileName);

  boost::iostreams::filtering_ostream fout;
  fout.push(boost::iostreams::gzip_compressor());
  fout.push(rawFile);

  const auto& basis = ssl_.getBasis();
  size_t nStates = basis.size();
  size_t steps   = duration() / interval + 1;

  std::vector<size_t> indices;

  for(size_t i = 0; i < nStates; ++i)
  {
    if(std::find(momNames.begin(), momNames.end(), basis[i]->getName()) != momNames.end())
      indices.push_back(i);
  }

  fout << "Generation";

  for(auto const& nm : momNames)
    fout << "," << nm;

  fout << "\n";

  auto matVar = engine_->toEigenMatrixVariant();
  auto vecWrap = MatrixEngine::VectorVariant(engine_->getVectorVariant()); // copy by value not to modify engine_

  visitMatrixAndVector(matVar, vecWrap, [&](auto const& M, auto& vWrap)
  {
    auto multiplyInPlace = [&](auto* vecPtr)
    {
      for(size_t k = 0; k < interval; ++k)
        *vecPtr = M * (*vecPtr);
    };

    for(size_t step = 0; step < steps; ++step)
    {
      long long generation = static_cast<long long>(startGen_) - static_cast<long long>(step) * interval;

      fout << generation;

      for(size_t idx : indices)
        fout << "," << static_cast<double>(vWrap[idx]);

      fout << "\n";

      if(step + 1 < steps)
        multiplyInPlace(&vWrap);
    }
  });

  fout.reset();
  rawFile.close();
}

//------------------------------------------------------------------------------
// Prints recursion formulas Δ[state] = … for each basis function except “I”
//------------------------------------------------------------------------------
void Epoch::printRecursions(std::ostream& stream)
{
  stream << '\n';
  const auto& basis = ssl_.getBasis();

  for(size_t i = 0; i < basis.size(); ++i)
  {
    auto const& state = basis[i];
    if(state->getName() == "I")
    continue;

    Eigen::Index pos = state->getPosition();
    stream << "\u0394[" << state->getName() << "] = ";

    bool firstTerm = true;

    for(auto const& op : operators_)
    {
      auto const& params = op->getParameters();

      for(size_t k = 0; k < params.size(); ++k)
      {
        auto const& param = params[k];
        double pval = param.getValue();

        if(pval == 0.0)
        continue;

        auto eigenVar = op->getMatrix(k).toEigenMatrixVariant();

        std::visit(overloaded{[&](auto const& M)
        {
          using MatT = std::decay_t<decltype(M)>;
          using Scalar = typename MatT::Scalar;

          Scalar scale = static_cast<Scalar>(pval);
          Eigen::Index ncols = M.cols();

          for(Eigen::Index col = 0; col < ncols; ++col)
          {
            Scalar coeff = M.coeff(pos, col) / scale;

            if(coeff == Scalar(0))
              continue;

            if(!firstTerm && coeff > Scalar(0))
              stream << '+';

            stream << std::fixed << std::setprecision(3) << coeff << '*' << param.getName() << '*' << basis[size_t(col)]->getName() << ' ';

            firstTerm = false;
          }
        }
        }, eigenVar);
      }
    }

    stream << '\n';
  }
}

//------------------------------------------------------------------------------
// Dumps the (summed) transition matrix to a file as row col value triples
//------------------------------------------------------------------------------
void Epoch::printTransitionMat(const std::string& fileName) const
{
  auto eigenVar = engine_->toEigenMatrixVariant();

  std::ofstream out(fileName, std::ios::binary);
  if(!out)
    throw bpp::Exception("printTransitionMat: cannot open " + fileName);

  // preserve caller’s formatting
  auto oldFlags = out.flags();
  auto oldPrec = out.precision();

  out << std::scientific << std::setprecision(16);

  std::visit(overloaded{[&](auto const& M)
  {
    using MatT   = std::decay_t<decltype(M)>;

    out << M.rows() << ' ' << M.cols() << ' ' << M.nonZeros() << '\n';

    for(Eigen::Index k = 0; k < M.outerSize(); ++k)
    {
      for(typename MatT::InnerIterator it(M, k); it; ++it)
        out << it.row() << ' ' << it.col() << ' ' << it.value() << '\n';
    }
  }
  }, eigenVar);

  // restores formatting
  out.flags(oldFlags);
  out.precision(oldPrec);
}

//------------------------------------------------------------------------------
// Computes steady state by largest eigenpair and update engine & moments
//------------------------------------------------------------------------------
void Epoch::computeEigenSteadyState()
{
  EigenResult res = findLeadingEigenpair();

  auto matVar = engine_->getMatrixVariant();
  MatrixEngine::VectorVariant vecWrap;

  std::visit(overloaded{[&](auto const& M)
  {
    using Scalar = typename std::decay_t<decltype(M)>::Scalar;

    // allocates an Eigen::Matrix<Scalar,Dynamic,1>
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> v(M.cols());

    for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(M.cols()); ++i)
      v(i) = static_cast<Scalar>(res.vector(i));

    vecWrap = Vector<Scalar>(std::move(v));
  }
  }, matVar);

  engine_->setVector(std::move(vecWrap));
  updateMoments(engine_->getVectorVariant());
}

//------------------------------------------------------------------------------
// Pseudo‐steady state via power‐method for discrete‐time model
//------------------------------------------------------------------------------
void Epoch::computePseudoSteadyStateDiscrete(double tol)
{
  bool converged = false;

  auto matVar = engine_->getMatrixVariant();
  MatrixEngine::VectorVariant vecWrap = engine_->getVectorVariant();

  #ifdef DEBUG
  auto demangle = [](const std::type_info& ti)
  {
    int status;
    char* demangled = abi::__cxa_demangle(ti.name(), nullptr, nullptr, &status);
    std::string result = (status == 0 && demangled) ? demangled : ti.name();
    free(demangled);
    return result;
  };

  std::visit([&](auto const& x)
  {
    std::cout << "matVar holds: " << demangle(typeid(x)) << "\n";
    std::cout << "matVar type hash: " << typeid(x).hash_code() << "\n";
    std::cout << "matVar object address: " << static_cast<const void*>(&x) << "\n";
  }, matVar);

  std::visit([&](auto const& x)
  {
    std::cout << "vecWrap holds: " << demangle(typeid(x)) << "\n";
    std::cout << "vecWrap type hash: " << typeid(x).hash_code() << "\n";
    std::cout << "vecWrap object address: " << static_cast<const void*>(&x) << "\n";
  }, vecWrap);
  #endif

  visitMatrixAndVector(matVar, vecWrap,[&](auto const& A, auto& dstWrap)
  {
    //A.print(std::cout);
    using Scalar = typename std::decay_t<decltype(A)>::Scalar;

    size_t pop = ssl_.getPopIndices()[0];
    double mu = getParameterValue("u_" + bpp::TextTools::toString(pop));
    double s = getParameterValue("s_" + bpp::TextTools::toString(pop));
    size_t twoN = static_cast<size_t>(1.0 / getParameterValue("1/2N_" + bpp::TextTools::toString(pop)));

    double hr = twoN * mu;
    double hl = (std::abs(s) < 1e-14) ? twoN * mu : (2. * twoN * mu * std::exp(2. * twoN * s) / (std::exp(2. * twoN * s) - 1.0) - mu/s);

    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> y(A.cols()); // inits y (Eigen::Vec)

    double f = 1.0; // scaling factor to mimic the factors of 1-2p
    for(Eigen::Index i = 0; i < y.size(); ++i)
    {
      const auto& prefix = ssl_.getBasis()[size_t(i)]->getPrefix();

      if(i > 0 && prefix == ssl_.getBasis()[size_t(i - 1)]->getPrefix()) // same kind of moment, applies heuristic decay
        f *= 0.925;

      else // resets for next moment
        f = 1.0;

      if(prefix == "Hl")
        y(i) = Scalar(hl * f);

      else if(prefix == "Hr")
        y(i) = Scalar(hr * f);

      else if(prefix == "pi2")
        y(i) = Scalar(hr * hl * f * 1.3);

      else if(prefix == "I")
        y(i) = Scalar(1.0);

      else if(prefix == "DD")
        y(i) = Scalar(hr * hl * f * 5e-1);

      else // if(prefix == "Dr")
        y(i) = Scalar(hr * hl * f * 3e-1);
    }

    // burn‐in
    size_t burnSteps = twoN;
    for(size_t b = 0; b < burnSteps; ++b)
    {
      Vector<Scalar> yWrap(std::move(y)); // wraps thin Eigen vector
      Vector<Scalar> result = A * yWrap;
      y = result.eigen();
    }

    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> prev = y;

    size_t maxIter = 20 * twoN;
    for(size_t iter = 0; iter < maxIter; ++iter)
    {
      //y = A * y;
      Vector<Scalar> yWrap(std::move(y)); // wraps thin Eigen vector
      Vector<Scalar> result = A * yWrap;
      y = result.eigen();

      if(iter % 10000 == 0)
      {
        double maxRel = 0.0;

        for(Eigen::Index i = 0; i < y.size(); ++i)
        {
          double pv = static_cast<double>(prev(i));
          double cv = static_cast<double>(y(i));
          maxRel = std::max(maxRel, std::abs(cv - pv) / std::max(1.0, std::abs(pv)));
        }

        if(maxRel <= tol)
        {
          std::cout << "\nPseudo steady-state converged after " << iter << " iterations, " << name_ << "\n";
          converged = true;
          break;
        }
      }

      prev = y;
    }

    // only write back into the wrapper if we converged
    if(converged)
      dstWrap = Vector<Scalar>(std::move(y));

  });

  if(!converged)
  {
    std::cerr << "\nPseudo steady-state did not converge. Falling back to eigen steady-state.\n";
    computeEigenSteadyState();
    return;
  }

  engine_->setVector(std::move(vecWrap));
  updateMoments(engine_->getVectorVariant());
}

//------------------------------------------------------------------------------
// Typed helper: does all the work in the Scalar domain
// Returns true if converged, writes steady‐state vector into outY
//------------------------------------------------------------------------------
template<typename Scalar>
bool Epoch::computePseudoSSContinuousTyped(
    const Matrix<Scalar>& A,
    double burnInTime,
    double dt,
    double tol,
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& outY) const
{
  const Eigen::Index n = A.cols();
  Eigen::Matrix<Scalar, Eigen::Dynamic, 1> y(n);

  // rough initial guess
  size_t pop = ssl_.getPopIndices()[0];
  double mu = getParameterValue("u_" + bpp::TextTools::toString(pop));
  double s = getParameterValue("s_" + bpp::TextTools::toString(pop));
  size_t twoN = static_cast<size_t>(1. / getParameterValue("1/2N_" + bpp::TextTools::toString(pop)));

  double hr = twoN * mu;
  double hl = (std::abs(s) < 1e-14) ? twoN * mu : (2. * twoN * mu * std::exp(2. * twoN * s) / (std::exp(2. * twoN * s) - 1.) - mu / s);

  double f = 1.0; // scaling factor to mimic the factors of 1-2p
  for(Eigen::Index i = 0; i < y.size(); ++i)
  {
    const auto& prefix = ssl_.getBasis()[size_t(i)]->getPrefix();

    if(i > 0 && prefix == ssl_.getBasis()[size_t(i - 1)]->getPrefix()) // same kind of moment, applies heuristic decay
      f *= 0.925;

    else // resets for next moment
      f = 1.0;

    if(prefix == "Hl")
      y(i) = Scalar(hl * f);

    else if(prefix == "Hr")
      y(i) = Scalar(hr * f);

    else if(prefix == "pi2")
      y(i) = Scalar(hr * hl * f * 1.3);

    else if(prefix == "I")
      y(i) = Scalar(1.0);

    else if(prefix == "DD")
      y(i) = Scalar(hr * hl * f * 5e-1);

    else // if(prefix == "Dr")
      y(i) = Scalar(hr * hl * f * 3e-1);
  }

  size_t burnSteps = static_cast<size_t>(burnInTime / dt);

  Vector<Scalar> Vwrap;

  for(size_t b = 0; b < burnSteps; ++b)
  {
    Vwrap.eigen() = y;

    if constexpr(std::is_same_v<Scalar, double>)
      y = integrateDoubleCN(A, Vwrap, dt, dt).eigen();

    else
      y = integrateMpfrCN(A, Vwrap, dt, dt).eigen();
  }

  Eigen::Matrix<Scalar, Eigen::Dynamic, 1> prev = y;
  size_t maxSteps = static_cast<size_t>(10 * burnInTime / dt);
  bool converged = false;

  for(size_t iter = 0; iter < maxSteps; ++iter)
  {
    Vwrap.eigen() = y;

    if constexpr(std::is_same_v<Scalar, double>)
      y = integrateDoubleCN(A, Vwrap, dt, dt).eigen();

    else
      y = integrateMpfrCN(A, Vwrap, dt, dt).eigen();

    if(iter % 10000 == 0)
    {
      double maxRel = 0.0;

      for(Eigen::Index i = 0; i < y.size(); ++i)
      {
        double pv = static_cast<double>(prev(i));
        double cv = static_cast<double>(y(i));
        maxRel = std::max(maxRel, std::abs(cv - pv) / std::max(1.0, std::abs(pv)));
      }

      if(maxRel <= tol)
      {
        std::cout << "\nPseudo steady-state converged after " << iter << " iterations, " << name_ << "\n";
        converged = true;
        break;
      }
    }

    prev = y;
  }

  if(!converged)
    return false;

  outY = std::move(y);

  return true;
}

//------------------------------------------------------------------------------
// Public dispatcher: unpack the variant and call the typed helper
//------------------------------------------------------------------------------
void Epoch::computePseudoSteadyStateContinuous(double burnInTime, double dt, double tol)
{
  auto matVar = engine_->getMatrixVariant();
  bool converged = false;

  MatrixEngine::VectorVariant vecWrap;  // will hold the wrapped steady‐state

  std::visit(overloaded{[&](auto const& A)
  {
    using Scalar = typename std::decay_t<decltype(A)>::Scalar;

    // Allocate an Eigen‐vector for the output
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> y;
    converged = computePseudoSSContinuousTyped(A, burnInTime, dt, tol, y);

    if(converged)
      vecWrap = Vector<Scalar>(std::move(y));
  }
  }, matVar);

  if(!converged)
  {
    std::cerr << "Crank–Nicolson steady-state did not converge. Falling back to eigen-based steady-state.\n";
    computeEigenSteadyState();
    return;
  }

  engine_->setVector(std::move(vecWrap));
  updateMoments(engine_->getVectorVariant());
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

//------------------------------------------------------------------------------
// Single‐precision vs. multi‐precision integrator dispatch
//------------------------------------------------------------------------------
MatrixEngine::VectorVariant
Epoch::integrate(double dt, double totalTime) const
{
  if(dt <= 0.0 || totalTime < 0.0)
    throw bpp::Exception("Epoch::integrate: dt must be >0 and totalTime >=0");

  const auto& matVar = engine_->getMatrixVariant();
  const auto& vecVar = engine_->getVectorVariant();
  MatrixEngine::VectorVariant result;

  visitMatrixAndVector(matVar, vecVar, [&](auto const& M, auto const& V)
  {
    using Scalar = typename std::decay_t<decltype(M)>::Scalar;

    if constexpr (std::is_same_v<Scalar, double>)
      result = integrateDoubleCN(M, V, dt, totalTime);

    else
      result = integrateMpfrCN(M, V, dt, totalTime);
  });

  return result;
}

//------------------------------------------------------------------------------
// Adaptive‐step integrator dispatch for double vs. mpreal
//------------------------------------------------------------------------------
MatrixEngine::VectorVariant
Epoch::integrateAdaptive(double dt,
                         double totalTime,
                         double tol,
                         double dtMin,
                         double dtMax) const
{
  if(dtMin <= 0.0 || dtMax <= 0.0 || dtMin > dtMax) // TODO improve condition checking
    throw bpp::Exception("Epoch::integrateAdaptive: invalid dtMin/dtMax");

  if(tol < 0.0)
    throw bpp::Exception("Epoch::integrateAdaptive: tol must be >=0");

  const auto& matVar = engine_->getMatrixVariant();
  const auto& vecVar = engine_->getVectorVariant();
  MatrixEngine::VectorVariant result;

  visitMatrixAndVector(matVar, vecVar, [&](auto const& M, auto const& V)
  {
    using Scalar = typename std::decay_t<decltype(M)>::Scalar;

    if constexpr(std::is_same_v<Scalar, double>)
      result = integrateAdaptiveDoubleCN(M, V, dt, totalTime, tol, dtMin, dtMax);

    else
      result = integrateAdaptiveMpfrCN(M, V, dt, totalTime, tol, dtMin, dtMax);

  });

  return result;
}

//------------------------------------------------------------------------------
// Build the full transition matrix and steady‐state engine
//------------------------------------------------------------------------------
void Epoch::init_()
{
  if(operators_.empty())
    throw bpp::Exception("Epoch::init_() called with no operators.");

  if(!engine_)
    throw bpp::Exception("Epoch::init_() called with null engine_.");

  #ifdef DEBUG
  auto demangle = [](const std::type_info& ti)
  {
    int status;
    char* demangled = abi::__cxa_demangle(ti.name(), nullptr, nullptr, &status);
    std::string result = (status == 0 && demangled) ? demangled : ti.name();
    free(demangled);
    return result;
  };

  operators_.front()->getParameters().printParameters(std::cout);
  #endif

  MatrixEngine::MatrixVariant accWrap = operators_.front()->getTransitionMatrixVariant();

  for(size_t i = 1; i < operators_.size(); ++i)
  {
    MatrixEngine::MatrixVariant nextWrap = operators_[i]->getTransitionMatrixVariant();

    #ifdef DEBUG
    operators_[i]->getParameters().printParameters(std::cout);

    std::visit([&](auto const& x)
    {
      std::cout << "accWrap holds: " << demangle(typeid(x)) << "\n";
      std::cout << "accWrap type hash: " << typeid(x).hash_code() << "\n";
      std::cout << "accWrap object address: " << static_cast<const void*>(&x) << "\n";
      //x.print(std::cout);
    }, accWrap);

    std::visit([&](auto const& x)
    {
      std::cout << "nextWrap holds: " << demangle(typeid(x)) << "\n";
      std::cout << "nextWrap type hash: " << typeid(x).hash_code() << "\n";
      std::cout << "nextWrap object address: " << static_cast<const void*>(&x) << "\n";
      //x.print(std::cout);
    }, nextWrap);
    #endif

    visitSameType(accWrap, nextWrap, [&](auto& A, auto const& B)
    {
      A += B;  // A and B are the same Matrix<T> type
    });

    #ifdef DEBUG
    std::visit([&](auto const& x)
    {
      std::cout << "POST-AGG accWrap holds: " << demangle(typeid(x)) << "\n";
      std::cout << "POST-AGG accWrap type hash: " << typeid(x).hash_code() << "\n";
    }, accWrap);
    #endif
  }

  engine_->setMatrix(std::move(accWrap));
  engine_->addIdentityInPlace();
  engine_->pruneInPlace();
  engine_->compressInPlace();
}

