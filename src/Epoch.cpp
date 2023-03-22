/*
 * Authors: Gustavo V. Barroso
 * Created: 31/08/2022
 * Last modified: 20/03/2023
 *
 */

#include <ios>

#include "Log.hpp"
#include "Epoch.hpp"

void Epoch::fireParameterChanged(const bpp::ParameterList& params)
{
  if(matchParametersValues(params))
    updateOperators_(params);

  Eigen::SparseMatrix<double> mat = operators_[0]->getTransitionMatrix(); // init mat

  for(size_t i = 1; i < operators_.size(); ++i)
    mat = mat * operators_[i]->getTransitionMatrix();

  transitionMatrix_ = mat; // converts from sparse to dense format
}

void Epoch::computeExpectedSumStats(Eigen::VectorXd& y)
{
  y = transitionMatrix_.pow(duration()) * y;
}

std::vector<size_t> Epoch::fetchSelectedPopIds()
{
  std::vector<size_t> ret(0);
  ret.reserve(pops_.size());

  for(auto it = std::begin(pops_); it != std::end(pops_); ++it)
  {
    if(it->second->hasSelection())
      ret.emplace_back(it->first);
  }

  return ret;
}

void Epoch::transferStatistics(Eigen::VectorXd& y) // y comes from previous Epoch
{
  Eigen::VectorXd tmp(ssl_.getCompressedBasis().size()); // y and tmp have potentially different sizes
  tmp.setZero();

  // for each Moment in *this Epoch, we assign its value from its parental Moment from the previous Epoch
  for(int i = 0; i < tmp.size(); ++i)
  {
    int idx = ssl_.getCompressedBasis()[i]->getParent()->getPosition();
    tmp(i) = y(idx);
  }

  y = tmp; // swap
}

void Epoch::updateMoments(const Eigen::VectorXd& y)
{
  assert(y.size() == static_cast<int>(ssl_.getCompressedBasis().size()));

  for(int i = 0; i < y.size(); ++i)
    ssl_.getCompressedBasis()[i]->setValue(y(i));
}

void Epoch::printRecursions(std::ostream& stream)
{
  for(size_t i = 0; i < ssl_.getCompressedBasis().size(); ++i)
  {
    if(ssl_.getCompressedBasis()[i]->getName() != "I")
    {
      int pos = static_cast<int>(ssl_.getCompressedBasis()[i]->getPosition()); // row in delta matrix
      stream << "\u0394[" << ssl_.getCompressedBasis()[i]->getName() << "] = ";

      for(size_t j = 0; j < operators_.size(); ++j)
      {
        for(size_t k = 0; k < operators_[j]->getParameters().size(); ++k)
        {
          bpp::Parameter param = operators_[j]->getParameters()[k];
          std::string name = param.getName();
          auto mat = operators_[j]->getMatrix(k); // hard copy delta matrix
          mat = mat / param.getValue(); // convert back to coefficients

          for(int l = 0; l < mat.cols(); ++l)
          {
            if(mat.coeffRef(pos, l) != 0)
            {
              if(mat.coeffRef(pos, l) > 0)
                stream << "+";

              stream << bpp::TextTools::toString(mat.coeffRef(pos, l)) + "*" + name + "*" + ssl_.getCompressedBasis()[l]->getName() + " ";
            }
          }
        }
      }

      stream << "\n";
    }
  }
}

void Epoch::computeSteadyState_()
{
  #ifdef VERBOSE
  Log logger;
  logger.openFile(getName() + "_matrices.txt");
  Eigen::SparseMatrix<double> test(ssl_.getNumCompressedStats(), ssl_.getNumCompressedStats());
  test.setIdentity();

  for(size_t i = 0; i < operators_.size(); ++i)
  {
    Eigen::SparseMatrix<double> tmp(ssl_.getNumCompressedStats(), ssl_.getNumCompressedStats());
    tmp.setZero();
    for(size_t j = 0; j < operators_[i]->getMatrices().size(); ++j)
    {
      logger.getLogFile() << "\n\nsum of entries (delta matrix " << j << ") = " << std::setprecision(1e-12) << std::scientific << operators_[i]->getMatrices()[j].sum() << "\n";
      tmp += operators_[i]->getMatrices()[j];
      bpp::ParameterList pl;
      pl.addParameter(operators_[i]->getParameters()[j]);
      pl.printParameters(logger.getLogFile());
      logger.getLogFile() << std::setprecision(1e-12) << std::scientific << operators_[i]->getMatrices()[j] << "\n";
    }

    operators_[i]->getParameters().printParameters(logger.getLogFile());
    logger.getLogFile() << "\n\nsum of entries (operator combined delta matrix) = " << std::setprecision(1e-12) << std::scientific << tmp.sum() << "\n";
    logger.getLogFile() << std::setprecision(1e-12) << tmp << "\n";
    logger.getLogFile() << "operator transition matrix:\n";
    logger.getLogFile() << operators_[i]->getTransitionMatrix() << "\n\n";

    logger.getLogFile() << "accumulated transition matrix:\n";
    test = operators_[i]->getTransitionMatrix() * test;
    logger.getLogFile() << std::setprecision(1e-12) << test << "\n";
  }
  #endif

  Eigen::SparseMatrix<double> mat = operators_[0]->getTransitionMatrix(); // init mat

  for(size_t i = 1; i < operators_.size(); ++i)
    mat = operators_[i]->getTransitionMatrix() * mat;

  transitionMatrix_ = mat; // converts from sparse to dense format
  Eigen::EigenSolver<Eigen::MatrixXd> es(transitionMatrix_);

  int idx = 0;
  for(int i = 0; i < es.eigenvalues().size(); ++i)
  {
    // finding the leading eigenvalue (should be == 1., but searching for equality is problematic due to precision)
    if(es.eigenvalues().real()(i) > es.eigenvalues().real()(idx))
      idx = i;
  }

  steadYstate_ = es.eigenvectors().col(idx).real();
  steadYstate_ /= steadYstate_(ssl_.findCompressedIndex(ssl_.getDummyIndexUncompressed())); // I moment embodies scaling constant used by Eigen

  updateMoments(steadYstate_);
}
