/*
 * Authors: Gustavo V. Barroso
 * Created: 31/08/2022
 * Last modified: 10/09/2025
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

void Epoch::transferStatistics(std::unique_ptr<VectorInterface>& y) // y comes from previous Epoch
{
  // y and tmp have potentially different sizes due to number of Populations and/or Order(1-2p)
  auto tmp = y->cloneWithSize(ssl_.getBasis().size()); // maintains the derived class of y

  // for each Moment in *this Epoch, we assign its value from its parental Moment from the previous Epoch
  // this follows the ancestry patterns of moments according to population history, see Model::linkMoments_()
  for(size_t i = 0; i < tmp->size(); ++i)
  {
    size_t parentPos = ssl_.getBasis()[i]->getParent()->getPosition();
    setValue(tmp.get(), i, getValue(y.get(), parentPos));
  }

  y = std::move(tmp);
}

void Epoch::updateMoments(std::unique_ptr<VectorInterface>& y)
{
  assert(y.size() == static_cast<int>(ssl_.getBasis().size()));

  for(int i = 0; i < y.size(); ++i)
    ssl_.getBasis()[i]->setValue(y(i).toDouble()); // WARNING must adopt type polymorhpism in SSL as well?
                                                   // NOTE do I ever extract the value from ssl_-> moments_ / basis_ for actual computation or is it there just for organization / printing?
}

void Epoch::printMoments(std::ostream& stream)
{
  std::vector<std::shared_ptr<Moment>> tmp = getSslib().getBasis();

  for(auto& m : tmp)
    stream << std::setprecision(24) << m->getName() << " = " << m->getValue() << "\n";
}

// prints expectations of Hl and Hr over time (precision not terribly importantly here)
void Epoch::printHetMomentsIntermediate(std::unique_ptr<VectorInterface>& y, const std::string& modelName, size_t interval)
{
  transferStatistics(y); // since different Epochs may use different Order

  std::string fileName = modelName + "_" + name_ + "_hets_time.txt";
  std::ofstream fout(fileName);

  std::vector<std::shared_ptr<Moment>> tmp = getSslib().getBasis();
  size_t steps = duration() / interval + 1; // prints every interval generations

  for(size_t i = 0; i < steps; ++i)
  {
    for(size_t j = 0; j < tmp.size(); ++j)
    {
      if(tmp[j]->getName() == "Hr_0_0" || tmp[j]->getName() == "Hl_0_0")
        fout << std::setprecision(16) << tmp[j]->getName() << " = " << y[j] << " " << startGen_ - i * interval << "\n"; // WARNING y[i]
    }

    if(i < steps - 1) { // not to advance further than needed, important when there are > 2 Epochs
      for(size_t k = 0; k < interval; ++k)
        y = transitionMatrix_ * y;
    }
  }

  fout.close();
}

void Epoch::printRecursions(std::ostream& stream)
{
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
          bpp::Parameter param = operators_[j]->getParameters()[k];
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
  transitionMatrix_ ->print(fileName);
}

void Epoch::computeEigenSteadyState()
{
  testSteadyState();
  init_();

  // converting to dense format to perform eigen decomposition
  Eigen::Matrix<mpfr::mpreal, Eigen::Dynamic,  Eigen::Dynamic> denseTransMat = transitionMatrix_;
  Eigen::EigenSolver<Eigen::Matrix<mpfr::mpreal, Eigen::Dynamic, Eigen::Dynamic>> es(denseTransMat);

  int idx = 0;
  for(int i = 0; i < es.eigenvalues().size(); ++i)
  {
    // finding the leading eigenvalue (== 1.,but not searching for equality due to precision)
    if(es.eigenvalues().real()(i) > es.eigenvalues().real()(idx))
      idx = i;
  }

  if(es.eigenvalues().real()(idx) > 1. + 1e-5) // NOTE
  {
    double cond = fetchConditionNumber();
    std::cout << "\nCondition Number of transition matrix = " << cond << "\n";
    std::cout << "\nLeading eigenvalue of full transition matrix = " << es.eigenvalues().real()(idx) << "\n";
    throw bpp::Exception("Epoch::Leading Eigenvalue > 1! Consider using a smaller order of 1-2p factors.\n");
  }

  // I moment embodies scaling constant used by Eigen
  steadYstate_ = es.eigenvectors().col(idx).real();
  steadYstate_ /= steadYstate_(ssl_.findCompressedIndex(ssl_.getMoment("I")));

  updateMoments(steadYstate_);
}

void Epoch::computePseudoSteadyState()
{
  testSteadyState();
  init_();

  if(dynamic_cast<MatrixDouble*>(transitionMatrix_.get()))
    steadYstate_ = std::make_unique<VectorDouble>(transitionMatrix_.size());

  else if(dynamic_cast<MatrixMPReal*>(transitionMatrix_.get()))
    steadYstate_ = std::make_unique<VectorMPReal>(transitionMatrix_.size());

  else
    throw bpp::Exception("AbstractOperator::Mis-cast transition matrix!");

  // a very rough guess for starting values to help w/ convergence
  size_t pop = ssl_.getPopIndices()[0];
  double mu = getParameterValue("u_" + bpp::TextTools::toString(pop));
  double s = getParameterValue("s_" + bpp::TextTools::toString(pop));
  int twoN = static_cast<int>(1 / getParameterValue("1/2N_" + bpp::TextTools::toString(pop)));
  double h = twoN * mu;

  for(int i = 0; i < y.size(); ++i)
  {
    const std::string& prefix = ssl_.getBasis()[i]->getPrefix();

    ifprefix == "Hl")
    {
      double gamma = static_cast<double>(twoN * s);

      if(gamma < -5.)
      {
        double p = mu / -s;
        transitionMatrix_(i) = p * (1-p);
      }

      else
      {
        transitionMatrix_(i) = h;
      }
    }

    else if(prefix == "Hr")
      transitionMatrix_(i) = h;

    else if(prefix == "pi2")
      transitionMatrix_(i) = h * h;

    else if(prefix == "I")
      transitionMatrix_(i) = 1.;

    else // DD and Dr
      transitionMatrix_(i) = h * h * 1e-2;
  }
  
  for(size_t j = 0; j < 10 * twoN; ++j)
    y = transitionMatrix_ * y;

  steadYstate_ = y;
  updateMoments(steadYstate_);
}

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

  // "split" original matrix in two
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
}

