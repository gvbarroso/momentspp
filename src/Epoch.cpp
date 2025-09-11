/*
 * Authors: Gustavo V. Barroso
 * Created: 31/08/2022
 * Last modified: 11/09/2025
 *
 */

#include <ios>

#include "Migration.hpp"
#include "Epoch.hpp"

namespace
{
  inline mpfr::mpreal getValue(const VectorInterface& vec, size_t index)
  {
    if(const auto* mpVec = dynamic_cast<const VectorMPReal*>(&vec))
      return mpVec->get(index);

    else
      return mpfr::mpreal(vec.get(index));
  }

  inline void setValue(VectorInterface& vec, size_t index, const mpfr::mpreal& value)
  {
    if(auto* mpVec = dynamic_cast<VectorMPReal*>(&vec))
      mpVec->set(index, value);

    else
      vec.set(index, value.toDouble());
  }

  inline void setValue(VectorInterface& vec, size_t index, double value)
  {
    vec.set(index, value);  // works for both types
  }
}

void Epoch::fireParameterChanged(const bpp::ParameterList& params)
{
  if(matchParametersValues(params))
  {
    updateOperators_(params);
    std::unique_ptr<MatrixInterface> sum = operators_[0]->getTransitionMatrix()->clone();

    for(size_t i = 1; i < operators_.size(); ++i)
      sum->addInPlace(operators_[i]->getTransitionMatrix().get());

    transitionMatrix_ = std::move(sum);
  }
}

void Epoch::computeExpectedSumStats(std::unique_ptr<VectorInterface>& y)
{
  for(size_t i = 0; i < duration(); ++i)
    y = transitionMatrix_->multiply(y.get());

  updateMoments(y);
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
void Epoch::transferStatistics(std::unique_ptr<VectorInterface>& y) // y comes from previous Epoch
{
  // y and tmp have potentially different sizes due to number of Populations and/or Order(1-2p)
  auto tmp = y->cloneWithSize(ssl_.getBasis().size()); // maintains the derived class of y

  // for each Moment in *this Epoch, we assign its value from its parental Moment from the previous Epoch
  // this follows the ancestry patterns of moments according to population history, see Model::linkMoments_()
  for(size_t i = 0; i < tmp->size(); ++i)
  {
    size_t parentPos = ssl_.getBasis()[i]->getParent()->getPosition();
    setValue(tmp.get(), i, getValue(y.get(), parentPos)); // inline function above
  }

  y = std::move(tmp);
}

void Epoch::updateMoments(const std::unique_ptr<VectorInterface>& y)
{
  assert(y.size() == static_cast<int>(ssl_.getBasis().size()));

  // double precision is fine as expections in ssl_-> moments_ / basis_ exist just for organization / printing
  for(int i = 0; i < y.size(); ++i)
    ssl_.getBasis()[i]->setValue(y->get(i));

}

void Epoch::printMoments(std::ostream& stream)
{
  std::vector<std::shared_ptr<Moment>> tmp = getSslib().getBasis();

  for(auto& m : tmp)
    stream << std::setprecision(24) << m->getName() << " = " << m->getValue() << "\n";
}

// prints expectations of Hl and Hr over time
void Epoch::printHetMomentsIntermediate(std::unique_ptr<VectorInterface>& y, const std::string& modelName, size_t interval)
{
  // NOTE method could be adapted to take moment names as input
  transferStatistics(y); // since different Epochs may use different Order

  std::string fileName = modelName + "_" + name_ + "_hets_time.txt";
  std::ofstream fout(fileName);

  const std::vector<std::shared_ptr<Moment>>& basis = getSslib().getBasis();
  size_t steps = duration() / interval + 1; // prints every interval generations

  for(size_t i = 0; i < steps; ++i)
  {
    for(size_t j = 0; j < basis.size(); ++j)
    {
      // precision not terribly importantly here
      if(basis[j]->getName() == "Hr_0_0" || basis[j]->getName() == "Hl_0_0")
        fout << basis[j]->getName() << " = " << y->get(j) << " " << startGen_ - i * interval << "\n";
    }

    if(i < steps - 1) { // not to advance further than needed, important when there are > 2 Epochs
      for(size_t k = 0; k < interval; ++k)
        y = transitionMatrix_->multiply(y.get());
    }
  }

  fout.close();
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

  // deducing type of steadYstate_ from the type of transitionMatrix_
  if(dynamic_cast<MatrixMPReal*>(transitionMatrix_.get()))
    steadYstate_ = std::make_unique<VectorMPReal>(res.vector);

  else if(dynamic_cast<MatrixDouble*>(transitionMatrix_.get()))
  {
    Eigen::VectorXd vecDouble(res.vector.size());
    for(Eigen::Index i = 0; i < res.vector.size(); ++i)
      vecDouble(i) = res.vector(i).toDouble();

    steadYstate_ = std::make_unique<VectorDouble>(vecDouble);
  }

  // I moment embodies scaling constant used by Eigen
  //steadYstate_ = es.eigenvectors().col(idx).real();
  //steadYstate_ /= steadYstate_(ssl_.findCompressedIndex(ssl_.getMoment("I")));
  updateMoments(steadYstate_);
}

// assumes discrete-time treatment is adequate
void Epoch::computePseudoSteadyState()
{
  if(dynamic_cast<MatrixDouble*>(transitionMatrix_.get()))
    steadYstate_ = std::make_unique<VectorDouble>(transitionMatrix_.size());

  else if(dynamic_cast<MatrixMPReal*>(transitionMatrix_.get()))
    steadYstate_ = std::make_unique<VectorMPReal>(transitionMatrix_.size());

  else
    throw bpp::Exception("Epoch::Mis-cast transition matrix!");

  size_t pop = ssl_.getPopIndices()[0];
  double mu = getParameterValue("u_" + bpp::TextTools::toString(pop));
  double s = getParameterValue("s_" + bpp::TextTools::toString(pop));
  size_t twoN = static_cast<size_t>(1 / getParameterValue("1/2N_" + bpp::TextTools::toString(pop)));

  // rough guesses for starting values to help w/ convergence
  double hr = twoN * mu;
  double hl = hr;

  if(static_cast<double>(twoN * s) < -5)
  {
    double p = mu / -s;
    hl = p * (1-p);
  }

  for(size_t i = 0; i < steadYstate_.size(); ++i)
  {
    const std::string& prefix = ssl_.getBasis()[i]->getPrefix();

    if(prefix == "Hl")
      steadYstate_->set(i, hl);

    else if(prefix == "Hr")
      steadYstate_->set(i, hr);

    else if(prefix == "pi2")
      steadYstate_->set(i, hr * hl);

    else if(prefix == "I")
      steadYstate_->set(i, 1.);

    else // DD and Dr
      steadYstate_->set(i, hr * hl * 1e-1);
  }
  
  double tol = 1e-3;
  auto prev = steadYstate_->clone();  // deep copy if needed

  auto notConverged = [&](size_t i)
  {
    double prevVal = prev->get(i);
    double currVal = steadYstate_->get(i);
    double relDiff = std::abs(currVal - prevVal) / std::max(1.0, std::abs(prevVal));
    return relDiff > tol;
  };

  while(true)
  {
    bool converged = true;
    for(size_t i = 0; i < steadYstate_->size(); ++i)
    {
      if(notConverged(i))
      {
        converged = false;
        break;
      }
    }

    if(converged)
      break;

    prev = steadYstate_->clone();  // update previous guess
    steadYstate_ = transitionMatrix_->multiply(steadYstate_.get());
  }

  updateMoments(steadYstate_);
}

// in models with gene-flow
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

std::unique_ptr<VectorInterface> Epoch::integrate(std::unique_ptr<VectorInterface> moms, double dt, size_t steps) const
{
  if(!transitionMatrix_ || !steadyState_) {
    throw bpp::Exception("Epoch not initialized!\n");

  auto I = transitionMatrix_->identity();
  auto halfDtM = transitionMatrix_->clone();

  // "splits" original matrix in two
  auto halfDtM->scale(dt / 2.0);
  auto A = I->add(*halfDtM->scale(-1.0));  // A = I - dt/2 * M
  auto B = I->add(*halfDtM);               // B = I + dt/2 * M

  std::unique_ptr<VectorInterface> y = moms->clone();

  for(size_t i = 0; i < steps; ++i)
  {
    auto rhs = B->multiply(*y);
    y = A->solve(*rhs);
  }

  return y;
}

void Epoch::init_()
{
  auto mat = operators_[0]->getTransitionMatrix(); // "delta" matrix

  for(size_t i = 1; i < operators_.size(); ++i) {
    mat->add(operators_[i]->getTransitionMatrix()); // summing rather than multiplying together
    mat.makeCompressed();
  }

  mat->add(operators_[0]->getIdentity()); // sums Identity to convert from "delta" to transition matrix
  mat.prune(0.0);  // removes converted zeros
  mat.makeCompressed();
  transitionMatrix_ = mat;

  //testSteadyState();
}

