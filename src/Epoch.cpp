/*
 * Authors: Gustavo V. Barroso
 * Created: 31/08/2022
 * Last modified: 15/09/2025
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
  if(matchParametersValues(params))
  {
    updateOperators_(params);
    engine_->setMatrix(operators_[0]->getTransitionMatrix());

    for(size_t i = 1; i < operators_.size(); ++i)
      engine_->addMatrixInPlace(operators_[i]->getTransitionMatrix());
  }
}

void Epoch::computeExpectedSumStatsDiscrete(const MatrixEngine::VectorVariantEigen& y)
{
  auto result = std::visit([&](const auto& vec) -> decltype(vec)
  {
    using Scalar = typename std::decay_t<decltype(vec)>::Scalar;

    const Eigen::SparseMatrix<Scalar> M = std::get<Eigen::SparseMatrix<Scalar>>(engine_->getRawMatrix());
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> y = vec;

    for(size_t i = 0; i < duration(); ++i)
      y = M * y;

    return y;
  }, y);

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
void Epoch::transferStatistics(MatrixEngine::VectorVariantEigen& y) const
{
  std::visit([&](auto& vecPrev)
  {
    using VectorType = std::decay_t<decltype(vecPrev)>;
    using Scalar = typename VectorType::Scalar;

    // y and tmp have potentially different sizes due to number of Populations and/or Order(1-2p)
    const size_t newSize = ssl_.getBasis().size();
    VectorType vecNew(newSize);

    // for each Moment in *this Epoch, we assign its value from its parental Moment from the previous Epoch
    // this follows the ancestry patterns of moments according to population history, see Model::linkMoments_()
    for(size_t i = 0; i < newSize; ++i)
    {
      size_t parentPos = ssl_.getBasis()[i]->getParent()->getPosition();
      vecNew(i) = vecPrev(parentPos);
    }

    vecPrev = std::move(vecNew);  // overwrite original vector in-place
  }, y);
}

void Epoch::updateMoments(const MatrixEngine::VectorVariantEigen& y)
{
  std::visit([&](const auto& vec)
  {
    using VectorType = std::decay_t<decltype(vec)>;
    assert(vec.size() == static_cast<int>(ssl_.getBasis().size()));

    for(int i = 0; i < vec.size(); ++i)
      ssl_.getBasis()[i]->setValue(static_cast<double>(vec(i)));

  }, y);
}

void Epoch::printMoments(std::ostream& stream)
{
  std::vector<std::shared_ptr<Moment>> tmp = getSslib().getBasis();

  for(auto& m : tmp)
    stream << std::setprecision(16) << m->getName() << " = " << m->getValue() << "\n";
}

void Epoch::printMomentsIntermediate(
  MatrixEngine::VectorVariantEigen& y,
  const std::string& modelName,
  size_t interval,
  const std::vector<std::string>& momNames)
{
  transferStatistics(y);  // adjust vector to match current Epoch's basis

  std::string fileName = modelName + "_" + name_ + "_moments.csv.gz";
  std::ofstream rawFile(fileName, std::ios_base::out | std::ios_base::binary);
  boost::iostreams::filtering_ostream fout;
  fout.push(boost::iostreams::gzip_compressor());
  fout.push(rawFile);

  const auto& basis = getSslib().getBasis();
  const size_t steps = duration() / interval + 1; // prints every interval generations

  // column indices for chosen moments
  std::vector<size_t> indices;
  for(size_t i = 0; i < basis.size(); ++i)
  {
    const std::string& name = basis[i]->getName();
    if(std::find(momNames.begin(), momNames.end(), name) != momNames.end())
      chosen.push_back(i);
  }

  fout << "Generation";
  for(const auto& name : momNames)
    fout << "," << name;
  fout << "\n";

  std::visit([&](auto& vec)
  {
    using VectorType = std::decay_t<decltype(vec)>;
    using Scalar = typename VectorType::Scalar;

    for(size_t i = 0; i < steps; ++i)
    {
      fout << (startGen_ - i * interval);  // gen column

      for(size_t idx : chosen)
        fout << "," << static_cast<double>(vec(idx));  // mom values

      fout << "\n";

      if(i < steps - 1) // not to advance further than needed, important when there are > 2 Epochs
      {
        for (size_t k = 0; k < interval; ++k)
        {
          const auto& M = std::get<Eigen::SparseMatrix<Scalar>>(engine_->getRawMatrix());
          vec = M * vec;
        }
      }
    }
  }, y);

  fout.reset();  // flush and close gzip stream
  rawFile.close();
}

void Epoch::printRecursions(std::ostream& stream)
{
  // this method only prints first-order coefficients,
  // even if transitionMatrix_ is obtained by multiplying operators
  stream << "\n";

  for(size_t i = 0; i < ssl_.getBasis().size(); ++i)
  {
    if(ssl_.getBasis()[i]->getName() != "I")
    {
      int pos = static_cast<int>(ssl_.getBasis()[i]->getPosition()); // row in delta matrix
      stream << "\u0394[" << ssl_.getBasis()[i]->getName() << "] = ";

      for(size_t j = 0; j < operators_.size(); ++j) // admixture coefficients are more complex
      {
        for(size_t k = 0; k < operators_[j]->getParameters().size(); ++k)
        {
          const bpp::Parameter& param = operators_[j]->getParameters()[k];
          const std::string& name = param.getName();

          auto mat = operators_[j]->getMatrix(k); // hard copy delta matrix

          if(param.getValue() != 0.)
            mat = mat / param.getValue(); // convert back to coefficients

          for(int l = 0; l < mat.cols(); ++l)
          {
            if(mat.coeffRef(pos, l) != 0)
            {
              if(mat.coeffRef(pos, l) > 0)
                stream << "+";

              stream << std::setprecision(3) << mat.coeffRef(pos, l) << "*" + name + "*" + ssl_.getBasis()[l]->getName() + " ";
            }
          }
        }
      }

      stream << "\n";
    }
  }
}

void Epoch::printTransitionMat(const std::string& fileName) const
{
  transitionMatrix_->print(fileName);
}

void Epoch::computeEigenSteadyState()
{
  EigenResult res = findLeadingEigenpair();

  if(res.value.toDouble() > 1. + 1e-14)
  {
    double cond = fetchConditionNumber();
    std::cout << "\nCondition Number of transition matrix = " << cond << "\n";
    std::cout << "\nLeading eigenvalue of full transition matrix = " << res.value.toDouble() << "\n";
    throw bpp::Exception("Epoch::Leading Eigenvalue > 1! Consider using a smaller order of 1-2p factors.\n");
  }

  // Deduce scalar type from engine's matrix and construct steady state vector
  std::visit([&](const auto& mat) {
    using Scalar = typename std::decay_t<decltype(mat)>::Scalar;

    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> converted(res.vector.size());

    for(Eigen::Index i = 0; i < res.vector.size(); ++i)
      converted(i) = static_cast<Scalar>(res.vector(i));

    engine_->setVector(MatrixEngine::VectorVariant{Vector<Scalar>(converted)});
  }, engine_->getMatrixVariant());

  updateMoments(engine_->getVectorVariant());
}

// assumes discrete-time treatment is adequate
void Epoch::computePseudoSteadyState(double tol = 1e-6)
{
  bool converged = false;

  std::visit([&](const auto& mat) {
    using Scalar = typename std::decay_t<decltype(mat)>::Scalar;

    Vector<Scalar> y(mat.mat_.rows());

    size_t pop = ssl_.getPopIndices()[0];
    double mu = getParameterValue("u_" + bpp::TextTools::toString(pop));
    double s = getParameterValue("s_" + bpp::TextTools::toString(pop));
    size_t twoN = static_cast<size_t>(1 / getParameterValue("1/2N_" + bpp::TextTools::toString(pop)));

    size_t maxIterations = 20 * twoN;

    double hr = twoN * mu;
    double hl = hr;

    if(static_cast<double>(twoN * s) < -5)
    {
      double p = mu / -s;
      hl = p * (1 - p);
    }

    for(size_t i = 0; i < y->size(); ++i)
    {
      const std::string& prefix = ssl_.getBasis()[i]->getPrefix();

      if (prefix == "Hl")
        y.set(i, Scalar(hl));
      else if (prefix == "Hr")
        y.set(i, Scalar(hr));
      else if (prefix == "pi2")
        y.set(i, Scalar(hr * hl));
      else if (prefix == "I")
        y.set(i, Scalar(1.0));
      else
        y.set(i, Scalar(hr * hl * 1e-1));
    }

    auto prev = y->clone();

    auto notConverged = [&](size_t i) {
      double prevVal = prev->get(i).toDouble();
      double currVal = y->get(i).toDouble();
      double relDiff = std::abs(currVal - prevVal) / std::max(1.0, std::abs(prevVal));
      return relDiff > tol;
    };

    for(size_t iter = 0; iter < maxIterations; ++iter)
    {
      bool allConverged = true;
      for(size_t i = 0; i < y->size(); ++i)
      {
        if(notConverged(i))
        {
          allConverged = false;
          break;
        }
      }

      if(allConverged)
      {
        converged = true;
        break;
      }

      prev = y->clone();
      y = mat.multiply(y);
    }

    if(converged)
     engine_->setVector(y);
  }, engine_->getMatrixVariant());

  if(!converged)
  {
    std::cerr << "Epoch::Pseudo steady state did not converge. Falling back to eigen-based steady state.\n";
    computeEigenSteadyState();
    return;
  }

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

template<typename Scalar>
Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Epoch::integrateTyped(const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& moms) const
{
  using MatrixType = Eigen::SparseMatrix<Scalar>;
  using VectorType = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

  const MatrixType M = std::get<MatrixType>(engine_->getRawMatrix());

  MatrixType I(M.rows(), M.cols());
  I.setIdentity();

  MatrixType halfDtM = dt_ / 2.0 * M;
  MatrixType A = I - halfDtM;
  MatrixType B = I + halfDtM;

  VectorType y = moms;

  Eigen::SparseLU<MatrixType> solver;
  solver.compute(A);
  if(solver.info() != Eigen::Success)
    throw bpp::Exception("Epoch::Matrix decomposition failed during integration.");

  for(size_t i = 0; i < steps_; ++i)
  {
    VectorType rhs = B * y;
    y = solver.solve(rhs);
    if (solver.info() != Eigen::Success)
      throw bpp::Exception("Epoch::Linear solve failed during integration.");
  }

  return y;
}

MatrixEngine::VectorVariantEigen Epoch::integrate(const MatrixEngine::VectorVariantEigen& moms) const
{
  return std::visit([&](const auto& vec) -> MatrixEngine::VectorVariantEigen
  {
    using Scalar = typename std::decay_t<decltype(vec)>::Scalar;
    return integrateTyped<Scalar>(vec);
  }, moms);
}

template<typename Scalar>
Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Epoch::integrateAdaptiveTyped(
  const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& moms,
  double dtMin = 1e-6,
  double dtMax = 1.0) const
{
  using MatrixType = Eigen::SparseMatrix<Scalar>;
  using VectorType = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

  MatrixType M = std::get<MatrixType>(engine_->getRawMatrix());

  MatrixType I(M.rows(), M.cols());
  I.setIdentity();

  VectorType y = moms;
  double t = 0.0;
  double dt = dt_;

  while(t < totalTime_)
  {
    if(t + dt > totalTime_)
      dt = totalTime_ - t;

    MatrixType halfDtM = dt / 2.0 * M;
    MatrixType A = I - halfDtM;
    MatrixType B = I + halfDtM;

    Eigen::SparseLU<MatrixType> solver;
    solver.compute(A);

    if(solver.info() != Eigen::Success)
      throw bpp::Exception("Epoch::Matrix decomposition failed during adaptive integration.");

    VectorType y_full = solver.solve(B * y);

    // Two half steps
    double halfDt = dt / 2.0;
    MatrixType halfDtM_half = halfDt * M;
    MatrixType A_half = I - halfDtM_half;
    MatrixType B_half = I + halfDtM_half;

    Eigen::SparseLU<MatrixType> solver_half;
    solver_half.compute(A_half);

    if(solver_half.info() != Eigen::Success)
      throw bpp::Exception("Epoch::Matrix decomposition failed for half step.");

    VectorType y_half1 = solver_half.solve(B_half * y);
    VectorType y_half2 = solver_half.solve(B_half * y_half1);

    double error = (y_full - y_half2).norm() / y_half2.norm();

    if(error < tolerance_)
    {
      y = y_half2;
      t += dt;
      dt = std::min(dt * 1.5, dtMax);
    }

    else
      dt = std::max(dt / 2.0, dtMin);
  }

  return y;
}

MatrixEngine::VectorVariantEigen Epoch::integrateAdaptive(const MatrixEngine::VectorVariantEigen& moms) const
{
  return std::visit([&](const auto& vec) -> MatrixEngine::VectorVariantEigen
  {
    using Scalar = typename std::decay_t<decltype(vec)>::Scalar;
    return integrateAdaptiveTyped<Scalar>(vec);
  }, moms);
}

void Epoch::init_()
{
  if(operators_.empty())
    throw bpp::Exception("Epoch::init_() called with no operators.");

  engine_->setMatrix(operators_[0]->getTransitionMatrix().getMatrixVariant());
  for(size_t i = 1; i < operators_.size(); ++i)
    engine_->addMatrixInPlace(operators_[i]->getTransitionMatrix().getMatrixVariant());

  engine_->addIdentityInPlace(); // sums Identity to convert from "delta" to transition matrix
  engine_->pruneInPlace();
  engine_->compressInPlace();
  // testSteadyState();
}


