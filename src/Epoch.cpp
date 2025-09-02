/*
 * Authors: Gustavo V. Barroso
 * Created: 31/08/2022
 * Last modified: 18/06/2024
 *
 */

#include <ios>

#include "Migration.hpp"
#include "Epoch.hpp"

void Epoch::fireParameterChanged(const bpp::ParameterList& params)
{
  if(matchParametersValues(params))
  {
    updateOperators_(params);

    Eigen::SparseMatrix<long double> mat = operators_[0]->getTransitionMatrix(); // init

    for(size_t i = 1; i < operators_.size(); ++i)
      mat = mat * operators_[i]->getTransitionMatrix();

    transitionMatrix_ = mat; // converts from sparse to dense format
  }
}

void Epoch::computeExpectedSumStats(Eigen::Matrix<long double, Eigen::Dynamic, 1>& y)
{
  // heavy linear algebra, uses Eigen multi-threading
  y = transitionMatrix_.pow(duration()) * y;
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

void Epoch::transferStatistics(Eigen::Matrix<long double, Eigen::Dynamic, 1>& y) // y comes from previous Epoch
{
  Eigen::Matrix<long double, Eigen::Dynamic, 1> tmp(ssl_.getBasis().size()); // y and tmp have potentially different sizes
  tmp.setZero();

  // for each Moment in *this Epoch, we assign its value from its parental Moment from the previous Epoch
  for(int i = 0; i < tmp.size(); ++i)
    tmp(i) = y(ssl_.getBasis()[i]->getParent()->getPosition());

  y = tmp;
}

void Epoch::updateMoments(const Eigen::Matrix<long double, Eigen::Dynamic, 1>& y)
{
  assert(y.size() == static_cast<int>(ssl_.getBasis().size()));

  for(int i = 0; i < y.size(); ++i)
    ssl_.getBasis()[i]->setValue(y(i));
}

void Epoch::printMoments(std::ostream& stream)
{
  std::vector<std::shared_ptr<Moment>> tmp = getSslib().getBasis();

  for(auto& m : tmp)
    stream << std::setprecision(24) << m->getName() << " = " << m->getValue() << "\n";
}

// prints expectations of Hl and Hr over time
void Epoch::printHetMomentsIntermediate(Eigen::Matrix<long double, Eigen::Dynamic, 1>& y, const std::string& modelName, size_t interval)
{
  transferStatistics(y); // since different Epochs may use different Order

  std::string fileName = modelName + "_" + name_ + "_hets_time.txt";
  std::ofstream fout(fileName);

  std::vector<std::shared_ptr<Moment>> tmp = getSslib().getBasis();
  size_t numTimeSteps = duration() / interval + 1; // prints every interval generations

  for(size_t i = 0; i < numTimeSteps; ++i)
  {
    for(size_t j = 0; j < tmp.size(); ++j)
    {
      if(tmp[j]->getName() == "Hr_0_0" || tmp[j]->getName() == "Hl_0_0")
        fout << std::setprecision(24) << tmp[j]->getName() << " = " << y[j] << " " << startGen_ - i * interval << "\n";
    }

    if(i < numTimeSteps - 1) // not to advance further than needed, important when there are > 2 Epochs
      y = transitionMatrix_.pow(interval) * y;
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
          std::string name = param.getName();

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
  std::ofstream matFile;
  matFile.open(fileName);

  for(int i = 0; i < transitionMatrix_.rows(); ++i)
  {
    for(int j = 0; j < transitionMatrix_.cols(); ++j)
    {
      matFile << transitionMatrix_.coeffRef(i, j);

      if(j < transitionMatrix_.cols() - 1)
        matFile << ",";
    }

    matFile  << "\n";
  }

  matFile.close();
}

void Epoch::computeEigenSteadyState()
{
  testSteadyState();
  init_();
  Eigen::EigenSolver<Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic>> es(transitionMatrix_);

  int idx = 0;
  for(int i = 0; i < es.eigenvalues().size(); ++i)
  {
    // finding the leading eigenvalue (== 1.,but not searching for equality due to precision)
    if(es.eigenvalues().real()(i) > es.eigenvalues().real()(idx))
      idx = i;
  }

  if(es.eigenvalues().real()(idx) > 1. + 1e-5)
  {
    double cond = fetchConditionNumber();
    std::cout << "\nCondition Number of transition matrix = " << cond << "\n";
    throw bpp::Exception("Epoch::Leading Eigenvalue > 1! Consider using a smaller order of 1-2p factors.\n");
  }

  // I moment embodies scaling constant used by Eigen
  steadYstate_ = es.eigenvectors().col(idx).real();
  steadYstate_ /= steadYstate_(ssl_.findCompressedIndex(ssl_.getMoment("I")));

  updateMoments(steadYstate_);
}

void Epoch::computePseudoSteadyState() // for speed
{
  testSteadyState();
  init_();

  Eigen::Matrix<long double, Eigen::Dynamic, 1> y(transitionMatrix_.rows());

  // a very rough guess for starting values to help w/ convergence
  size_t p = ssl_.getPopIndices()[0];
  long double h = getParameterValue("u_" + bpp::TextTools::toString(p)) / getParameterValue("1/2N_" + bpp::TextTools::toString(p));

  for(int i = 0; i < y.size(); ++i)
  {
    if(ssl_.getBasis()[i]->getPrefix() == "Hl" || ssl_.getBasis()[i]->getPrefix() == "Hr")
      y(i) = h;

    else if(ssl_.getBasis()[i]->getPrefix() == "pi2")
      y(i) = h * h * 1e-1;

    else if(ssl_.getBasis()[i]->getPrefix() == "I")
      y(i) = 1.;

    else
      y(i) = h * 1e-4;
  }

  steadYstate_ = transitionMatrix_.pow(1e+6) * y; // in practice 1e+6 gens. is good enough?
  updateMoments(steadYstate_);
}

void Epoch::testSteadyState()
{
  /*if(pops_.size() > 1)
  {
    for(size_t i = 0; i < operators_.size(); ++i)
    {
      auto tmp = std::dynamic_pointer_cast<Migration>(operators_[i]);

      if(tmp != nullptr)
        tmp->testFlow();
    }
  }*/
}

void Epoch::init_()
{
  Eigen::SparseMatrix<long double> mat = operators_[0]->getTransitionMatrix();

  for(size_t i = 1; i < operators_.size(); ++i)
    mat = mat * operators_[i]->getTransitionMatrix();

  transitionMatrix_ = mat; // converts from sparse to dense format
}

